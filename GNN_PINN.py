# -*- coding: utf-8 -*-
"""
3-Stage Flow:
(1) Non-eq data => Pretrain => MSE only
(2) eq data => separate saturP model => MSE only
(3) eq data => use predicted saturP => PINN Finetune
"""
   
import os
import numpy as np
import pandas as pd
import rdkit
import rdkit.Chem
from rdkit.Chem import AllChem
import tensorflow as tf
from tensorflow.keras.callbacks import EarlyStopping
import matplotlib.pyplot as plt

# CPU mode (optional)
os.environ["CUDA_VISIBLE_DEVICES"] = "0"


###############################################################################
# Hyperparameters
###############################################################################
EPOCHS_PRETRAIN = 1000
EPOCHS_SATUR    = 500
EPOCHS_FINETUNE = 1000

LR_PRETRAIN     = 1e-4
LR_SATUR        = 1e-4
LR_FINETUNE     = 1e-4

LAMBDA_PINN_PRETRAIN = 0.0   # Pretrain => no PINN
LAMBDA_PINN_FINETUNE = 0.5   # Finetune => PINN on

BATCH_SIZE = 16

early_stop = EarlyStopping(
    monitor='val_loss',
    patience=20,
    min_delta=0.0001,
    mode='min',
    restore_best_weights=True,
)

###############################################################################
# 1. CSV Load
###############################################################################
csv_path = "dataset.csv"
df_all = pd.read_csv(csv_path)
df_all = df_all.sample(frac=1, random_state=0).reset_index(drop=True)

# eq / non-eq 구분
df_eq  = df_all[(df_all["maskHv"] + df_all["maskSv"] + df_all["maskHl"] + df_all["maskSl"]) == 4]
df_non = df_all[(df_all["maskHv"] + df_all["maskSv"] + df_all["maskHl"] + df_all["maskSl"]) < 4]

def split_df(df, ratio=(0.8,0.1,0.1), seed=0):
    df = df.sample(frac=1, random_state=seed).reset_index(drop=True)
    n   = df.shape[0]
    r1, r2, r3 = ratio
    idx1= int(r1*n)
    idx2= int((r1+r2)*n)
    train_df= df.iloc[:idx1]
    val_df=   df.iloc[idx1:idx2]
    test_df=  df.iloc[idx2:]
    return train_df,val_df,test_df

df_non_train, df_non_val, df_non_test = split_df(df_non,(0.8,0.1,0.1), seed=0)
df_eq_train,  df_eq_val,  df_eq_test  = split_df(df_eq,  (0.8,0.1,0.1), seed=0)

###############################################################################
# 2. Normalization
###############################################################################
vars_to_norm= ["T","P","Hv","Hl","Sv","Sl","Tc","Pc","af"]

df_train_for_meanstd = pd.concat([df_non_train, df_eq_train], ignore_index=True)

mean_dict, std_dict= {}, {}
for c in vars_to_norm:
    mean_dict[c] = df_train_for_meanstd[c].mean()
    std_dict[c]  = df_train_for_meanstd[c].std() + 1e-12

def apply_zscore(df_sub):
    df_sub = df_sub.copy()
    for c in vars_to_norm:
        df_sub[c] = (df_sub[c] - mean_dict[c]) / std_dict[c]
    return df_sub

df_non_train_norm = apply_zscore(df_non_train)
df_non_val_norm   = apply_zscore(df_non_val)
df_non_test_norm  = apply_zscore(df_non_test)

df_eq_train_norm  = apply_zscore(df_eq_train)
df_eq_val_norm    = apply_zscore(df_eq_val)
df_eq_test_norm   = apply_zscore(df_eq_test)

###############################################################################
# 3. GNN/Utils
###############################################################################
my_elements = {
    6:  {"symbol":"C",  "value1":1.05, "value2":3.851, "value3":5.431, "value4":11.714, "value5":1.1346},
    8:  {"symbol":"O",  "value1":0.6,  "value2":3.5,   "value3":8.714, "value4":17.136, "value5":0.8597},
    1:  {"symbol":"H",  "value1":0.44, "value2":2.886, "value3":4.528, "value4":13.89,  "value5":0.8271},
    7:  {"symbol":"N",  "value1":0.69, "value2":3.66,  "value3":6.688, "value4":13.244, "value5":0.977},
    9:  {"symbol":"F",  "value1":0.5,  "value2":3.364, "value3":6.416, "value4":22.262, "value5":0.7686},
    16: {"symbol":"S",  "value1":2.74, "value2":4.035, "value3":6.928, "value4":8.972,  "value5":1.2806},
    17: {"symbol":"Cl", "value1":2.27, "value2":3.947, "value3":5.821, "value4":14.546, "value5":1.1615},
    35: {"symbol":"Br", "value1":2.51, "value2":4.189, "value3":5.692, "value4":17.52,  "value5":1.2982},
    53: {"symbol":"I",  "value1":3.39, "value2":4.5,   "value3":5.431, "value4":11.44,  "value5":1.527},
    15: {"symbol":"P",  "value1":3.05, "value2":4.147, "value3":5.463, "value4":8,      "value5":1.4362},
    34: {"symbol":"Se", "value1":2.91, "value2":4.205, "value3":6.428, "value4":8.262,  "value5":1.3906},
    50: {"symbol":"Sn", "value1":5.67, "value2":4.392, "value3":3.987, "value4":6.248,  "value5":1.8389},
    14: {"symbol":"Si", "value1":4.02, "value2":4.295, "value3":4.168, "value4":6.974,  "value5":1.6474},
    5:  {"symbol":"B",  "value1":1.05, "value2":3.851, "value3":5.431, "value4":11.714, "value5":1.1346},
    33: {"symbol":"As", "value1":5.67, "value2":4.392, "value3":3.987, "value4":6.248,  "value5":1.8389}
}
lookup = list(my_elements.keys())
max_graph_size = 150

def gen_smiles2graph(smiles):
    m= rdkit.Chem.MolFromSmiles(smiles)
    m= rdkit.Chem.AddHs(m)
    N= len(list(m.GetAtoms()))
    nodes= np.zeros((N, len(my_elements),5), dtype=np.float32)

    for atom in m.GetAtoms():
        idx= atom.GetIdx()
        atomic_num= atom.GetAtomicNum()
        if atomic_num not in lookup:
            continue
        elem_idx= lookup.index(atomic_num)
        for j in range(1,6):
            nodes[idx, elem_idx, j-1] = my_elements[atomic_num][f"value{j}"]

    bond_order = {
        rdkit.Chem.rdchem.BondType.SINGLE: 1.0,
        rdkit.Chem.rdchem.BondType.DOUBLE: 2.0,
        rdkit.Chem.rdchem.BondType.TRIPLE: 3.0,
        rdkit.Chem.rdchem.BondType.AROMATIC:1.5,
    }
    adj = np.zeros((N,N), dtype=np.float32)
    for bond in m.GetBonds():
        u= bond.GetBeginAtomIdx()
        v= bond.GetEndAtomIdx()
        val= bond_order.get(bond.GetBondType(), 0.0)
        adj[u,v]= val
        adj[v,u]= val
    adj+= np.eye(N,dtype=np.float32)

    return nodes, adj

def pad_array(arr, max_size):
    N= arr.shape[0]
    if N<max_size:
        pad_shape= (max_size-N,)+arr.shape[1:]
        arr= np.concatenate([arr, np.zeros(pad_shape,dtype=arr.dtype)], axis=0)
    else:
        arr= arr[:max_size]
    return arr

def pad_matrix(mat, max_size):
    N= mat.shape[0]
    if N<max_size:
        pad_mat= np.zeros((max_size,max_size), dtype=mat.dtype)
        pad_mat[:N,:N] = mat
        return pad_mat
    else:
        return mat[:max_size,:max_size]
    
class GCNLayer(tf.keras.layers.Layer):
    def __init__(self,activation=None,reduce_dim=False,**kwargs):
        super().__init__(**kwargs)
        self.activation= tf.keras.activations.get(activation)
        self.reduce_dim= reduce_dim

    def build(self,input_shape):
        node_shape, adj_shape= input_shape
        self.w= self.add_weight(shape=(len(my_elements),len(my_elements)),
                                name="w",
                                initializer=tf.random_normal_initializer())
        if self.reduce_dim:
            self.dense = self.add_weight(shape=(5,16),  name="dense")
            self.dense2= self.add_weight(shape=(16,16), name="dense2")
            self.dense3= self.add_weight(shape=(16,16), name="dense3")
            self.dense4= self.add_weight(shape=(16,),   name="dense4")

    def call(self,inputs):
        nodes, adj= inputs
        if self.reduce_dim:
            nodes= tf.einsum("bijk,km,mn,nl,l->bij",
                             nodes, self.dense, self.dense2, self.dense3, self.dense4)
        degree= tf.reduce_sum(adj, axis=-1)+1e-7
        new_nodes= tf.einsum("bi,bij,bjk,km->bim",1/degree, adj, nodes, self.w)
        out= self.activation(new_nodes)
        return out, adj

class GRLayer(tf.keras.layers.Layer):
    def call(self,inputs):
        nodes, adj= inputs
        mean_= tf.reduce_mean(nodes,axis=1)
        max_= tf.reduce_max(nodes,axis=1)
        min_= tf.reduce_min(nodes,axis=1)
        var_= tf.math.reduce_variance(nodes,axis=1)
        sum_= tf.reduce_sum(nodes,axis=1)
        lse_= tf.reduce_logsumexp(nodes,axis=1)

        c1= tf.expand_dims(lse_,axis=-1)
        c2= tf.expand_dims(mean_,axis=-1)
        c3= tf.expand_dims(max_,axis=-1)
        c4= tf.expand_dims(min_,axis=-1)
        c5= tf.expand_dims(var_,axis=-1)
        c6= tf.expand_dims(sum_,axis=-1)

        cat= tf.concat([c1,c2,c3,c4,c5,c6],axis=-1)
        shape= tf.shape(cat)
        new_shape=(shape[0], shape[1]*shape[2])
        return tf.reshape(cat,new_shape)

ninput= tf.keras.Input((None,len(my_elements),5), name="ninput")
ainput= tf.keras.Input((None,None),               name="ainput")

T_in=  tf.keras.Input(shape=(1,), name="T_in")
P_in=  tf.keras.Input(shape=(1,), name="P_in")
tr_in= tf.keras.Input(shape=(1,), name="tr_in")
pr_in= tf.keras.Input(shape=(1,), name="pr_in")
af_in= tf.keras.Input(shape=(1,), name="af_in")
Tc_in= tf.keras.Input(shape=(1,), name="Tc_in")
Pc_in= tf.keras.Input(shape=(1,), name="Pc_in")

x, A= GCNLayer("relu", reduce_dim=True)([ninput, ainput])
x, A= GCNLayer("relu")([x,A])
x= GRLayer()([x,A])

combined= tf.keras.layers.Concatenate()([
    x, T_in,P_in, tr_in, pr_in, af_in, Tc_in, Pc_in
])
x= tf.keras.layers.Dense(128,"relu")(combined)
x= tf.keras.layers.Dense(128,"relu")(x)
x= tf.keras.layers.Dense(64,"relu")(x)
x= tf.keras.layers.Dense(64,"relu")(x)

outputs= tf.keras.layers.Dense(4,activation=None)(x)

multi_model= tf.keras.Model(
    inputs=[ninput, ainput, T_in,P_in,tr_in,pr_in,af_in,Tc_in,Pc_in],
    outputs=outputs
)

def undo_norm_np(x, var_name):
    return x * std_dict[var_name] + mean_dict[var_name]

###############################################################################
# 4. Saturation Dataset => eq => (graph, T, af) => P
###############################################################################
def gen_satur_dataset(df_sub):
    for idx, row in df_sub.iterrows():
        smiles = str(row["SMILES"])
        nodes, adj= gen_smiles2graph(smiles)
        nodes= pad_array(nodes, max_graph_size)
        adj=   pad_matrix(adj,   max_graph_size)

        T_val= row["T"]   
        af_val= row.get("af",0.0)
        p_val= row["P"]   

        model_in= (
            nodes,
            adj,
            np.array([T_val], dtype=np.float32),
            np.array([af_val],dtype=np.float32)
        )
        yield (model_in, np.array([p_val],dtype=np.float32))

def make_satur_dataset(df_sub, batch_size=16):
    ds= tf.data.Dataset.from_generator(
        lambda: gen_satur_dataset(df_sub),
        output_signature=( 
            (
              tf.TensorSpec((max_graph_size,len(my_elements),5),tf.float32),
              tf.TensorSpec((max_graph_size,max_graph_size),tf.float32),
              tf.TensorSpec((1,),tf.float32),
              tf.TensorSpec((1,),tf.float32),
            ),
            tf.TensorSpec((1,),tf.float32)
        )
    )
    return ds.batch(batch_size)

###############################################################################
# 4. Non-eq Dataset => MSE
###############################################################################
def gen_non_dataset(df_sub):
    for idx, row in df_sub.iterrows():
        smiles = str(row["SMILES"])
        nodes, adj = gen_smiles2graph(smiles)
        nodes = pad_array(nodes, max_graph_size)
        adj   = pad_matrix(adj,   max_graph_size)

        T_val  = row["T"]    # (정규화된 T)
        P_val  = row["P"]    # (정규화된 P)
        tr_val = row["tr"]
        pr_val = row["pr"]
        af_val = row.get("af",0.0)
        Tc_val = row.get("Tc",1.0)
        Pc_val = row.get("Pc",1.0)

        Hv_ = row["Hv"]
        Sv_ = row["Sv"]
        Hl_ = row["Hl"]
        Sl_ = row["Sl"]

        mHv = row["maskHv"]
        mSv = row["maskSv"]
        mHl = row["maskHl"]
        mSl = row["maskSl"]

        # === 정규화 풀기(원래 P) ===
        P_orig = (row["P"] * std_dict["P"]) + mean_dict["P"]

        y_true = np.array([Hv_, Sv_, Hl_, Sl_], dtype=np.float32)
        y_mask = np.array([mHv,mSv,mHl,mSl],    dtype=np.float32)
        smiles_str = str(smiles)

        model_in = (
            nodes, adj,
            np.array([T_val], dtype=np.float32),
            np.array([P_val], dtype=np.float32),
            np.array([tr_val],dtype=np.float32),
            np.array([pr_val],dtype=np.float32),
            np.array([af_val],dtype=np.float32),
            np.array([Tc_val],dtype=np.float32),
            np.array([Pc_val],dtype=np.float32),
        )
        # label_in에 P_orig(1-D array 형태) 추가
        label_in = (
            y_true,
            y_mask,
            np.array([T_val],   dtype=np.float32),
            np.array([P_orig],  dtype=np.float32), 
            smiles_str,
        )
        yield (model_in, label_in)

def make_non_dataset(df_sub, batch_size=16):
    ds = tf.data.Dataset.from_generator(
        lambda: gen_non_dataset(df_sub),
        output_signature=( 
            (
                tf.TensorSpec((max_graph_size,len(my_elements),5), tf.float32),
                tf.TensorSpec((max_graph_size,max_graph_size),     tf.float32),
                tf.TensorSpec((1,), tf.float32),
                tf.TensorSpec((1,), tf.float32),
                tf.TensorSpec((1,), tf.float32),
                tf.TensorSpec((1,), tf.float32),
                tf.TensorSpec((1,), tf.float32),
                tf.TensorSpec((1,), tf.float32),
                tf.TensorSpec((1,), tf.float32),
            ),
            (
                tf.TensorSpec((4,), tf.float32),     # y_true
                tf.TensorSpec((4,), tf.float32),     # y_mask
                tf.TensorSpec((1,), tf.float32),     # T_val
                tf.TensorSpec((1,), tf.float32),     # P_orig
                tf.TensorSpec((),   tf.string),      # smiles_str
            )
        )
    )
    return ds.batch(batch_size)

###############################################################################
# 4. eqPINN Dataset => use "P_pred" in place of real P
###############################################################################
def gen_eqPINN_dataset(df_sub):
    """
    eq df 에서 'P_pred' 열은 PINN 입력용 (정규화된 satur P),
    'P' 열은 본래의 실측 P(정규화 전 or 후?)인데, 현재 df_sub는 z-score가 되어있으므로 
    원래 P = row["P"] * std_dict["P"] + mean_dict["P"] 로 복원
    """
    for idx, row in df_sub.iterrows():
        smiles = str(row["SMILES"])
        nodes, adj = gen_smiles2graph(smiles)
        nodes = pad_array(nodes, max_graph_size)
        adj   = pad_matrix(adj,   max_graph_size)

        T_val  = row["T"]         # z-score
        P_pred = row["P_pred"]    # z-score 상태 (saturation model 예측값)
        tr_val = row["tr"]
        pr_val = row["pr"]
        af_val = row.get("af",0.0)
        Tc_val = row.get("Tc",1.0)
        Pc_val = row.get("Pc",1.0)

        Hv_ = row["Hv"]
        Sv_ = row["Sv"]
        Hl_ = row["Hl"]
        Sl_ = row["Sl"]

        mHv = row["maskHv"]
        mSv = row["maskSv"]
        mHl = row["maskHl"]
        mSl = row["maskSl"]

        y_true = np.array([Hv_, Sv_, Hl_, Sl_], dtype=np.float32)
        y_mask = np.array([mHv,mSv,mHl,mSl],    dtype=np.float32)
        smiles_str= str(smiles)

        # === 원래 P ===
        P_orig = (row["P"] * std_dict["P"]) + mean_dict["P"]

        model_in = (
            nodes, adj,
            np.array([T_val],  dtype=np.float32),
            np.array([P_pred], dtype=np.float32),  # PINN은 satur P_pred(정규화) 입력
            np.array([tr_val], dtype=np.float32),
            np.array([pr_val], dtype=np.float32),
            np.array([af_val], dtype=np.float32),
            np.array([Tc_val], dtype=np.float32),
            np.array([Pc_val], dtype=np.float32),
        )
        # label_in에 P_orig 추가
        label_in = (
            y_true,
            y_mask,
            np.array([T_val],   dtype=np.float32),
            np.array([P_orig],  dtype=np.float32),  # 추가!
            smiles_str,
        )
        yield (model_in, label_in)

def make_eqPINN_dataset(df_sub, batch_size=16):
    ds = tf.data.Dataset.from_generator(
        lambda: gen_eqPINN_dataset(df_sub),
        output_signature=( 
            (
                tf.TensorSpec((max_graph_size,len(my_elements),5),tf.float32),
                tf.TensorSpec((max_graph_size,max_graph_size),tf.float32),
                tf.TensorSpec((1,),tf.float32),
                tf.TensorSpec((1,),tf.float32),
                tf.TensorSpec((1,),tf.float32),
                tf.TensorSpec((1,),tf.float32),
                tf.TensorSpec((1,),tf.float32),
                tf.TensorSpec((1,),tf.float32),
                tf.TensorSpec((1,),tf.float32),
            ),
            (
                tf.TensorSpec((4,),tf.float32),    # y_true
                tf.TensorSpec((4,),tf.float32),    # y_mask
                tf.TensorSpec((1,),tf.float32),    # T_val
                tf.TensorSpec((1,),tf.float32),    # P_orig
                tf.TensorSpec((),   tf.string),    # smiles_str
            )
        )
    )
    return ds.batch(batch_size)

###############################################################################
# 5. Saturation Model => (graph, T, af) => P
###############################################################################
sinput_n= tf.keras.Input((None,len(my_elements),5), name="ninput_sat")
sinput_a= tf.keras.Input((None,None),               name="ainput_sat")
sinput_T= tf.keras.Input(shape=(1,), name="T_sat")
sinput_af=tf.keras.Input(shape=(1,), name="af_sat")

sx, sA= GCNLayer("relu", reduce_dim=True)([sinput_n, sinput_a])
sx, sA= GCNLayer("relu")([sx, sA])
sx= GRLayer()([sx, sA])

sx_comb= tf.keras.layers.Concatenate()([sx, sinput_T, sinput_af])
sx= tf.keras.layers.Dense(128,"relu")(sx_comb)
sx= tf.keras.layers.Dense(64,"relu")(sx)
p_out= tf.keras.layers.Dense(1, activation=None)(sx)

model_satP= tf.keras.Model(
    inputs=[sinput_n, sinput_a, sinput_T, sinput_af],
    outputs= p_out
)

class SatPTrainer(tf.keras.Model):
    def __init__(self, sat_model):
        super().__init__()
        self.sat_model= sat_model
        self.loss_tracker= tf.keras.metrics.Mean(name="loss")

    def train_step(self, data):
        (model_in, p_true)= data
        with tf.GradientTape() as tape:
            p_pred= self.sat_model(model_in, training=True)
            loss_val= tf.reduce_mean(tf.square(p_pred- p_true))  # MSE
        grads= tape.gradient(loss_val, self.sat_model.trainable_variables)
        self.optimizer.apply_gradients(zip(grads, self.sat_model.trainable_variables))
        self.loss_tracker.update_state(loss_val)
        return {"loss": self.loss_tracker.result()}

    def test_step(self, data):
        (model_in, p_true)= data
        p_pred= self.sat_model(model_in, training=False)
        loss_val= tf.reduce_mean(tf.square(p_pred- p_true))
        self.loss_tracker.update_state(loss_val)
        return {"loss": self.loss_tracker.result()}

    @property
    def metrics(self):
        return [self.loss_tracker]

###############################################################################
# 6. PINN Model => multi_model
###############################################################################
ninput= tf.keras.Input((None,len(my_elements),5), name="ninput")
ainput= tf.keras.Input((None,None),               name="ainput")

T_in=  tf.keras.Input(shape=(1,), name="T_in")
P_in=  tf.keras.Input(shape=(1,), name="P_in")
tr_in= tf.keras.Input(shape=(1,), name="tr_in")
pr_in= tf.keras.Input(shape=(1,), name="pr_in")
af_in= tf.keras.Input(shape=(1,), name="af_in")
Tc_in= tf.keras.Input(shape=(1,), name="Tc_in")
Pc_in= tf.keras.Input(shape=(1,), name="Pc_in")

x, A= GCNLayer("relu", reduce_dim=True)([ninput, ainput])
x, A= GCNLayer("relu")([x,A])
x= GRLayer()([x,A])

combined= tf.keras.layers.Concatenate()([
    x, T_in,P_in,tr_in,pr_in,af_in,Tc_in,Pc_in
])
x= tf.keras.layers.Dense(128,"relu")(combined)
x= tf.keras.layers.Dense(128,"relu")(x)
x= tf.keras.layers.Dense(64,"relu")(x)
x= tf.keras.layers.Dense(64,"relu")(x)

outputs= tf.keras.layers.Dense(4, activation=None)(x)
multi_model= tf.keras.Model(
    inputs=[ninput, ainput, T_in,P_in,tr_in,pr_in,af_in,Tc_in,Pc_in],
    outputs=outputs
)

###############################################################################
# 7. Custom PINN Loss (BC와 HS 제거됨)
###############################################################################
def custom_pinn_loss(y_pred, y_true, mask, T_norm, P_norm, tr_norm, pr_norm, lambda_pinn):
    def undo_norm(tensor, var_name):
        return tensor * std_dict[var_name] + mean_dict[var_name]

    # y_pred => (batch,4) => (Hv, Sv, Hl, Sl) - 정규화 상태
    Hv_pred_n = y_pred[:,0]
    Sv_pred_n = y_pred[:,1]
    Hl_pred_n = y_pred[:,2]
    Sl_pred_n = y_pred[:,3]

    # y_true => (batch,4) => (Hv, Sv, Hl, Sl) - 정규화 상태
    Hv_true_n = y_true[:,0]
    Sv_true_n = y_true[:,1]
    Hl_true_n = y_true[:,2]
    Sl_true_n = y_true[:,3]

    # 역정규화
    Hv_pred = undo_norm(Hv_pred_n, "Hv")
    Sv_pred = undo_norm(Sv_pred_n, "Sv")
    Hl_pred = undo_norm(Hl_pred_n, "Hl")
    Sl_pred = undo_norm(Sl_pred_n, "Sl")

    Hv_true = undo_norm(Hv_true_n, "Hv")
    Sv_true = undo_norm(Sv_true_n, "Sv")
    Hl_true = undo_norm(Hl_true_n, "Hl")
    Sl_true = undo_norm(Sl_true_n, "Sl")

    # 필요 시 T의 역정규화
    T_real = undo_norm(tf.reshape(T_norm,[-1]), "T")
    P_real = undo_norm(tf.reshape(P_norm,[-1]), "P")

    # 1) MSE (마스크 적용)
    diff_Hv = (Hv_pred - Hv_true) * mask[:,0]
    diff_Sv = (Sv_pred - Sv_true) * mask[:,1]
    diff_Hl = (Hl_pred - Hl_true) * mask[:,2]
    diff_Sl = (Sl_pred - Sl_true) * mask[:,3]

    sq_err = diff_Hv**2 + diff_Sv**2 + diff_Hl**2 + diff_Sl**2
    denom = tf.reduce_sum(mask, axis=-1) + 1e-7
    sample_mse = sq_err / denom
    masked_mse = tf.reduce_mean(sample_mse)

    # 2) 평형 데이터에 한해 PDE 항목 계산
    sum_mask = tf.reduce_sum(mask, axis=-1)
    cond_eq = tf.equal(sum_mask, 4.0)

    Gv_ = Hv_pred - T_real * (Sv_pred/1000.0)
    Gl_ = Hl_pred - T_real * (Sl_pred/1000.0)
    pde_each = tf.square(Gv_ - Gl_)
    pde_each = tf.where(cond_eq, pde_each, tf.zeros_like(pde_each))
    pde_val = tf.reduce_mean(pde_each)

    # 3) 최종 Loss: MSE + lambda_pinn * PDE
    total_loss = masked_mse + lambda_pinn * pde_val

    return total_loss, masked_mse, pde_val

###############################################################################
# 8. Trainer classes
###############################################################################
class MultiHeadPINNTrainer(tf.keras.Model):
    def __init__(self, base_model, lambda_pinn):
        super().__init__()
        self.base_model= base_model
        self.loss_tracker= tf.keras.metrics.Mean(name="loss")
        self.lambda_pinn= lambda_pinn

    def train_step(self, data):
        (model_in, (y_true, y_mask, T_val, P_orig, smiles_str))= data
        with tf.GradientTape() as tape:
            y_pred= self.base_model(model_in, training=True)
            T_ = model_in[2]  # T_norm
            P_ = model_in[3]  # P_norm (in eq stage => P_pred)
            tr_= model_in[4]
            pr_= model_in[5]
            total_loss, mse_val, pde_val = custom_pinn_loss(
                y_pred, y_true, y_mask, T_, P_, tr_, pr_, self.lambda_pinn
            )
        grads= tape.gradient(total_loss, self.base_model.trainable_variables)
        self.optimizer.apply_gradients(zip(grads, self.base_model.trainable_variables))
        self.loss_tracker.update_state(total_loss)
        return {
            "loss": self.loss_tracker.result(),
            "mse":  mse_val,
            "pde":  pde_val,
        }

    def test_step(self, data):
        (model_in, (y_true,y_mask,T_val,P_orig,smiles_str))= data
        y_pred= self.base_model(model_in, training=False)
        T_= model_in[2]
        P_= model_in[3]
        tr_= model_in[4]
        pr_= model_in[5]

        total_loss, mse_val, pde_val = custom_pinn_loss(
            y_pred, y_true, y_mask, T_, P_, tr_, pr_, self.lambda_pinn
        )
        self.loss_tracker.update_state(total_loss)
        return {
            "loss": self.loss_tracker.result(),
            "mse":  mse_val,
            "pde":  pde_val,
        }

    @property
    def metrics(self):
        return [self.loss_tracker]

class SatPTrainer(tf.keras.Model):
    def __init__(self, sat_model):
        super().__init__()
        self.sat_model= sat_model
        self.loss_tracker= tf.keras.metrics.Mean(name="loss")

    def train_step(self, data):
        (model_in, p_true)= data
        with tf.GradientTape() as tape:
            p_pred= self.sat_model(model_in, training=True)
            loss_val= tf.reduce_mean(tf.square(p_pred- p_true))  # MSE
        grads= tape.gradient(loss_val, self.sat_model.trainable_variables)
        self.optimizer.apply_gradients(zip(grads, self.sat_model.trainable_variables))
        self.loss_tracker.update_state(loss_val)
        return {"loss": self.loss_tracker.result()}

    def test_step(self, data):
        (model_in, p_true)= data
        p_pred= self.sat_model(model_in, training=False)
        loss_val= tf.reduce_mean(tf.square(p_pred- p_true))
        self.loss_tracker.update_state(loss_val)
        return {"loss": self.loss_tracker.result()}

    @property
    def metrics(self):
        return [self.loss_tracker]

###############################################################################
# 9. Pretrain => non-eq => MSE only
###############################################################################
trainer_pre= MultiHeadPINNTrainer(base_model=multi_model, lambda_pinn=LAMBDA_PINN_PRETRAIN)
trainer_pre.compile(optimizer=tf.keras.optimizers.Adam(LR_PRETRAIN))

pretrain_train_ds= make_non_dataset(df_non_train_norm, BATCH_SIZE)
pretrain_val_ds=   make_non_dataset(df_non_val_norm,   BATCH_SIZE)

# dummy pass => create variables
dummy_nodes= np.zeros((1,max_graph_size,len(my_elements),5),dtype=np.float32)
dummy_adj=   np.zeros((1,max_graph_size,max_graph_size),dtype=np.float32)
dummy_T=     np.zeros((1,1),dtype=np.float32)
dummy_P=     np.zeros((1,1),dtype=np.float32)
dummy_tr=    np.zeros((1,1),dtype=np.float32)
dummy_pr=    np.zeros((1,1),dtype=np.float32)
dummy_af=    np.zeros((1,1),dtype=np.float32)
dummy_Tc=    np.zeros((1,1),dtype=np.float32)
dummy_Pc=    np.zeros((1,1),dtype=np.float32)

_ = multi_model([
    dummy_nodes, dummy_adj,
    dummy_T, dummy_P,
    dummy_tr, dummy_pr,
    dummy_af, dummy_Tc, dummy_Pc
])

print("Model variables created (Pretrain).")
if os.path.isdir("pretrain_tf"):
    print("Loading pretrain_tf ...")
    trainer_pre.load_weights("pretrain_tf")
else:
    print("No pretrain => train from scratch (non-eq MSE).")
    hist_pre= trainer_pre.fit(
        pretrain_train_ds,
        validation_data= pretrain_val_ds,
        epochs= EPOCHS_PRETRAIN,
        callbacks=[early_stop]
    )
    trainer_pre.save_weights("pretrain_tf")

###############################################################################
# 10. Saturation Model => eq => MSE
###############################################################################
sinput_n= tf.keras.Input((None,len(my_elements),5), name="ninput_sat")
sinput_a= tf.keras.Input((None,None),               name="ainput_sat")
sinput_T= tf.keras.Input(shape=(1,), name="T_sat")
sinput_af=tf.keras.Input(shape=(1,), name="af_sat")

sx, sA= GCNLayer("relu", reduce_dim=True)([sinput_n, sinput_a])
sx, sA= GCNLayer("relu")([sx, sA])
sx= GRLayer()([sx, sA])

sx_comb= tf.keras.layers.Concatenate()([sx, sinput_T, sinput_af])
sx= tf.keras.layers.Dense(128,"relu")(sx_comb)
sx= tf.keras.layers.Dense(64,"relu")(sx)
p_out= tf.keras.layers.Dense(1)(sx)

model_satP= tf.keras.Model(
    inputs=[sinput_n, sinput_a, sinput_T, sinput_af],
    outputs= p_out
)

trainer_sat= SatPTrainer(model_satP)
trainer_sat.compile(optimizer=tf.keras.optimizers.Adam(LR_SATUR))

sat_train_ds= make_satur_dataset(df_eq_train_norm, BATCH_SIZE)
sat_val_ds  = make_satur_dataset(df_eq_val_norm,   BATCH_SIZE)

if os.path.isdir("saturP_tf"):
    print("Loading saturP_tf ...")
    trainer_sat.load_weights("saturP_tf")
else:
    print("No saturP => train eq satur model(MSE).")
    hist_sat= trainer_sat.fit(
        sat_train_ds,
        validation_data= sat_val_ds,
        epochs= EPOCHS_SATUR,
        callbacks=[early_stop]
    )
    trainer_sat.save_weights("saturP_tf")

# satur test ds if needed
sat_test_ds= make_satur_dataset(df_eq_test_norm, BATCH_SIZE)
test_satur= trainer_sat.evaluate(sat_test_ds)
print(f"Satur Model Test Loss= {test_satur}")

###############################################################################
# 11. Predict saturP => eq df => 'P_pred'
###############################################################################
p_train= model_satP.predict(sat_train_ds).flatten()
df_eq_train_norm["P_pred"]= p_train

p_val= model_satP.predict(sat_val_ds).flatten()
df_eq_val_norm["P_pred"]= p_val

p_test= model_satP.predict(sat_test_ds).flatten()
df_eq_test_norm["P_pred"]= p_test

###############################################################################
# 12. Finetune => eq => PINN with P_pred
###############################################################################
finetune_train_ds= make_eqPINN_dataset(df_eq_train_norm, BATCH_SIZE)
finetune_val_ds=   make_eqPINN_dataset(df_eq_val_norm,   BATCH_SIZE)

trainer_fine= MultiHeadPINNTrainer(base_model=multi_model, lambda_pinn=LAMBDA_PINN_FINETUNE)
trainer_fine.compile(optimizer=tf.keras.optimizers.Adam(LR_FINETUNE))

_ = multi_model([
    dummy_nodes,dummy_adj,
    dummy_T, dummy_P,
    dummy_tr,dummy_pr,
    dummy_af,dummy_Tc,dummy_Pc
])
print("Variables created (Finetune).")

if os.path.isdir("finetune_tf"):
    print("Loading finetune_tf ...")
    trainer_fine.load_weights("finetune_tf")
else:
    print("No finetune => run eq PINN with saturP.")
    hist_fine= trainer_fine.fit(
        finetune_train_ds,
        validation_data= finetune_val_ds,
        epochs= EPOCHS_FINETUNE,
        callbacks=[early_stop]
    )
    trainer_fine.save_weights("finetune_tf")

###############################################################################
# 13. eq test => PINN with predicted P => Save CSV
###############################################################################
# 1) Make test dataset (equilibrium + P_pred)
finetune_test_ds = make_eqPINN_dataset(df_eq_test_norm, BATCH_SIZE)

# 2) Evaluate
test_result = trainer_fine.evaluate(finetune_test_ds, return_dict=True)
print("===== Finetune eq test (PINN with saturP) =====")
print(
  f"Loss={test_result['loss']:.4f}, "
  f"MSE={test_result['mse']:.4f}, PDE={test_result['pde']:.4f}"
)

# 다시 한 번 (같은) finetune_test_ds 만들어서, 예측 → CSV 저장
finetune_test_ds = make_eqPINN_dataset(df_eq_test_norm, BATCH_SIZE)

# --------------------------------------------
# 14. non-eq 테스트도 함께 예측 → 하나의 CSV
# --------------------------------------------
test_results = []  # eq/noneq 모두 담을 공용 리스트

###############################################################################
# (A) non-eq 예측
###############################################################################
for (model_inputs, (y_true, y_mask, T_val, P_orig, smiles_str)) in make_non_dataset(df_non_test_norm, BATCH_SIZE):
    y_pred = trainer_fine.base_model.predict(model_inputs, verbose=0)
    T_arr = model_inputs[2].numpy()  # (batch,1)
    P_arr = model_inputs[3].numpy()  # (batch,1)

    batch_size_ = y_pred.shape[0]

    for i in range(batch_size_):
        # -------------------------
        # 1) 예측값 (정규화 상태)
        # -------------------------
        Hv_pred_norm = y_pred[i, 0]
        Sv_pred_norm = y_pred[i, 1]
        Hl_pred_norm = y_pred[i, 2]
        Sl_pred_norm = y_pred[i, 3]

        # -------------------------
        # 2) 실제값 (정규화 상태)
        # -------------------------
        Hv_true_norm = y_true[i, 0]
        Sv_true_norm = y_true[i, 1]
        Hl_true_norm = y_true[i, 2]
        Sl_true_norm = y_true[i, 3]

        # 마스크
        maskHv_i = y_mask[i, 0]
        maskSv_i = y_mask[i, 1]
        maskHl_i = y_mask[i, 2]
        maskSl_i = y_mask[i, 3]

        # -------------------------
        # 3) T, P 역정규화
        # -------------------------
        T_norm = float(T_arr[i, 0])
        P_norm = float(P_arr[i, 0])

        T_real = T_norm * std_dict["T"] + mean_dict["T"]
        P_real = P_norm * std_dict["P"] + mean_dict["P"]

        # label에 들어있던 '원래 P'
        P_original = float(P_orig[i, 0])

        # -------------------------
        # 4) 예측값 역정규화 → float
        # -------------------------
        if float(maskHv_i) == 1.0:
            Hv_pred_real = float(Hv_pred_norm * std_dict["Hv"] + mean_dict["Hv"])
        else:
            Hv_pred_real = None
        
        if float(maskSv_i) == 1.0:
            Sv_pred_real = float(Sv_pred_norm * std_dict["Sv"] + mean_dict["Sv"])
        else:
            Sv_pred_real = None
        
        if float(maskHl_i) == 1.0:
            Hl_pred_real = float(Hl_pred_norm * std_dict["Hl"] + mean_dict["Hl"])
        else:
            Hl_pred_real = None
        
        if float(maskSl_i) == 1.0:
            Sl_pred_real = float(Sl_pred_norm * std_dict["Sl"] + mean_dict["Sl"])
        else:
            Sl_pred_real = None

        # -------------------------
        # 5) 실제값(ground truth)도 역정규화 → float
        # -------------------------
        if float(maskHv_i) == 1.0:
            Hv_true_real = float(Hv_true_norm * std_dict["Hv"] + mean_dict["Hv"])
        else:
            Hv_true_real = None
        
        if float(maskSv_i) == 1.0:
            Sv_true_real = float(Sv_true_norm * std_dict["Sv"] + mean_dict["Sv"])
        else:
            Sv_true_real = None
        
        if float(maskHl_i) == 1.0:
            Hl_true_real = float(Hl_true_norm * std_dict["Hl"] + mean_dict["Hl"])
        else:
            Hl_true_real = None
        
        if float(maskSl_i) == 1.0:
            Sl_true_real = float(Sl_true_norm * std_dict["Sl"] + mean_dict["Sl"])
        else:
            Sl_true_real = None

        # -------------------------
        # 6) 결과 저장
        # -------------------------
        smiles_py = smiles_str[i].numpy().decode("utf-8")

        test_results.append({
            "eq_type":       "non_eq",
            "SMILES":        smiles_py,
            "T_real":        T_real,
            "P_in_model":    P_real,
            "P_real_origin": P_original,
            "Hv_pred":       Hv_pred_real,
            "Sv_pred":       Sv_pred_real,
            "Hl_pred":       Hl_pred_real,
            "Sl_pred":       Sl_pred_real,
            "Hv_true":       Hv_true_real,
            "Sv_true":       Sv_true_real,
            "Hl_true":       Hl_true_real,
            "Sl_true":       Sl_true_real,
        })

###############################################################################
# (B) eq 예측
###############################################################################
for (model_inputs, (y_true, y_mask, T_val, P_orig, smiles_str)) in finetune_test_ds:
    y_pred = trainer_fine.base_model.predict(model_inputs, verbose=0)
    T_arr = model_inputs[2].numpy()
    P_arr = model_inputs[3].numpy()  # saturP_pred (정규화)

    batch_size_ = y_pred.shape[0]

    for i in range(batch_size_):
        # 1) 예측값
        Hv_pred_norm = y_pred[i, 0]
        Sv_pred_norm = y_pred[i, 1]
        Hl_pred_norm = y_pred[i, 2]
        Sl_pred_norm = y_pred[i, 3]

        # 2) 실제값
        Hv_true_norm = y_true[i, 0]
        Sv_true_norm = y_true[i, 1]
        Hl_true_norm = y_true[i, 2]
        Sl_true_norm = y_true[i, 3]

        maskHv_i = y_mask[i, 0]
        maskSv_i = y_mask[i, 1]
        maskHl_i = y_mask[i, 2]
        maskSl_i = y_mask[i, 3]

        # 3) T, saturP_pred 역정규화 → float
        T_norm = float(T_arr[i, 0])
        P_norm = float(P_arr[i, 0])
        T_real = T_norm * std_dict["T"] + mean_dict["T"]
        P_sat_pred = P_norm * std_dict["P"] + mean_dict["P"]

        # 실측 원래 P
        P_original = float(P_orig[i, 0])

        # 4) 예측값 역정규화
        if float(maskHv_i) == 1.0:
            Hv_pred_real = float(Hv_pred_norm * std_dict["Hv"] + mean_dict["Hv"])
        else:
            Hv_pred_real = None

        if float(maskSv_i) == 1.0:
            Sv_pred_real = float(Sv_pred_norm * std_dict["Sv"] + mean_dict["Sv"])
        else:
            Sv_pred_real = None

        if float(maskHl_i) == 1.0:
            Hl_pred_real = float(Hl_pred_norm * std_dict["Hl"] + mean_dict["Hl"])
        else:
            Hl_pred_real = None

        if float(maskSl_i) == 1.0:
            Sl_pred_real = float(Sl_pred_norm * std_dict["Sl"] + mean_dict["Sl"])
        else:
            Sl_pred_real = None

        # 5) 실제값 역정규화
        if float(maskHv_i) == 1.0:
            Hv_true_real = float(Hv_true_norm * std_dict["Hv"] + mean_dict["Hv"])
        else:
            Hv_true_real = None

        if float(maskSv_i) == 1.0:
            Sv_true_real = float(Sv_true_norm * std_dict["Sv"] + mean_dict["Sv"])
        else:
            Sv_true_real = None

        if float(maskHl_i) == 1.0:
            Hl_true_real = float(Hl_true_norm * std_dict["Hl"] + mean_dict["Hl"])
        else:
            Hl_true_real = None

        if float(maskSl_i) == 1.0:
            Sl_true_real = float(Sl_true_norm * std_dict["Sl"] + mean_dict["Sl"])
        else:
            Sl_true_real = None

        smiles_py = smiles_str[i].numpy().decode("utf-8")

        test_results.append({
            "eq_type":       "eq",
            "SMILES":        smiles_py,
            "T_real":        T_real,
            "P_in_model":    P_sat_pred,
            "P_real_origin": P_original,
            "Hv_pred":       Hv_pred_real,
            "Sv_pred":       Sv_pred_real,
            "Hl_pred":       Hl_pred_real,
            "Sl_pred":       Sl_pred_real,
            "Hv_true":       Hv_true_real,
            "Sv_true":       Sv_true_real,
            "Hl_true":       Hl_true_real,
            "Sl_true":       Sl_true_real,
        })

# --------------------------------------------
# (C) 최종 CSV 저장
# --------------------------------------------
df_results = pd.DataFrame(test_results)
df_results.to_csv("final_test_pred_all11111.csv", index=False)
print("Saved final_test_pred_all.csv.")
print(df_results.head(20))
print("===== Done Equilibrium Test & CSV Saving =====")
