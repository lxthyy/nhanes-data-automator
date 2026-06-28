#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NHANES 数据下载与处理工具 v2.03
================================
功能：一站式提取、清洗与验证 CDC NHANES 数据（1999-2024）

:author: 李鑫 (Li Xin)
:email: lxddzyx@126.com
:institution: 湖北医药学院附属太和医院
:license: MIT License
:version: 2.03
:repository: https://github.com/lxthyy/nhanes-downloader

Copyright (c) 2025 李鑫
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

数据来源：CDC NHANES (National Health and Nutrition Examination Survey)
"""

__author__ = "李鑫 (Li Xin)"
__email__ = "lxddzyx@126.com"
__license__ = "MIT"
__version__ = "2.03"

import os
import sys
import io
import json
import threading
import queue
import time
import urllib.request
import urllib.error
import ssl
import re
import struct
from datetime import datetime
from tkinter import (
    Tk, Frame, Label, Checkbutton, Button, Text, Scrollbar,
    filedialog, messagebox, simpledialog, ttk, BooleanVar, StringVar, Entry,
    OptionMenu, IntVar, Listbox, MULTIPLE, END, Canvas, PanedWindow
)
from tkinter import ttk

import pandas as pd
import numpy as np

# 0.1N 独立验证引擎
try:
    from validation_engine import ValidationEngine, gui_verify
    _VALIDATION_AVAILABLE = True
except ImportError:
    _VALIDATION_AVAILABLE = False

# ============================================================================
# 全局配置
# ============================================================================

NHANES_CYCLES = [
    ("1999-2000", 1999, "A"),
    ("2001-2002", 2001, "B"),
    ("2003-2004", 2003, "C"),
    ("2005-2006", 2005, "D"),
    ("2007-2008", 2007, "E"),
    ("2009-2010", 2009, "F"),
    ("2011-2012", 2011, "G"),
    ("2013-2014", 2013, "H"),
    ("2015-2016", 2015, "I"),
    ("2017-2018", 2017, "J"),
    ("2019-2020", 2019, "K"),
    ("2021-2022", 2021, "L"),
    ("2023-2024", 2023, "M"),
]

CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nhanes_cache")
os.makedirs(CACHE_DIR, exist_ok=True)

# ============================================================================
# 数据表别名 & 变量名兼容映射（不同周期表名/变量名有变化）
# ============================================================================
TABLE_ALIASES = {
    "THYROD": ["THYROD", "TST"],
    "TCHOL": ["TCHOL", "L13"],
    "GLU": ["GLU", "L10"],
    "INS": ["INS", "L10AM"],
    "BIOPRO": ["BIOPRO", "L16", "L40"],
    "HSCRP": ["HSCRP", "CRP"],
    "CBC": ["CBC", "L25"],
    "PBCD": ["PBCD", "L24", "L06"],
    "VID": ["VID", "VITB12", "VIC"],
    "UIO": ["UIO", "UMS"],
    "SEXHRM": ["SEXHRM", "TST"],
    "DR1TOT": ["DR1TOT", "DRXFMT"],
    "DR2TOT": ["DR2TOT", "DRXFMT2"],
    "RXQ_RX": ["RXQ_RX", "RXQ"],
}

VARIABLE_ALIASES = {
    # 甲状腺 - CDC 各周期变量名不同
    "LBXTSH": ["LBXTSH", "LBXTSH1", "LBDTSH1S", "LBDTSH"],
    "LBXT4": ["LBXT4", "LBXTT4", "LBXT4F", "LBDT4FSI", "LBDTT4S", "LBDTT4", "LBDT4FS"],
    "LBXFT4": ["LBXFT4", "LBXT4F", "LBDT4FS"],
    "LBXFT3": ["LBXFT3", "LBXT3F", "LBDFT3S", "LBDFT3"],
    "LBXTT3": ["LBXTT3", "LBDTT3SI", "LBDTT3S"],
    "LBXTG": ["LBXTG", "LBXTGN", "LBDTGNSI"],
    "LBXTPO": ["LBXTPO", "LBDTPOS"],
    "LBXATG": ["LBXATG", "LBXTGB"],
    # 血脂（2021+ LBD前缀）
    "LBXTC": ["LBXTC", "LBXSCH", "LBDTCS", "LBDTC"],
    "LBDHDL": ["LBDHDL", "LBXHDL", "LBDHDD", "LBDHDLS"],
    "LBDLDL": ["LBDLDL", "LBXLDL", "LBXSCLDL", "LBDLDLSI", "LBDLDLS"],
    "LBXTR": ["LBXTR", "LBXSTR", "LBDTRS", "LBDTR"],
    # 血糖（2021+ LBD前缀）
    "LBXGLU": ["LBXGLU", "LBXSGL", "LBDGLUS", "LBDGLU"],
    "LBXGH": ["LBXGH", "LBXGHB"],
    # 肌酐/尿酸（2021+ LBD前缀）
    "LBXSCR": ["LBXSCR", "LBXSCRSI", "LBDSCRS", "LBDSCR"],
    "LBXSUA": ["LBXSUA", "LBXSUASI", "LBDSUAS", "LBDSUA"],
    # 可替宁
    "LBXCOT": ["LBXCOT", "LBXCOTT"],
    # 肌酸激酶
    "LBXSCK": ["LBXSCK", "LBXSCKSI"],
    # 血常规
    "LBXWBC": ["LBXWBC", "LBXSWBC", "LBXWBCSI"],
    "LBXHGB": ["LBXHGB", "LBXSHGB"],
    "LBXPLT": ["LBXPLT", "LBXSPLT", "LBXPLTSI"],
    "LBXNEUT": ["LBXNEUT", "LBXNEPCT", "LBDNEUNO"],
    "LBXLY": ["LBXLY", "LBDLYMNO"],
    "LBXMONO": ["LBXMONO", "LBDMONO"],
    "LBXEOS": ["LBXEOS", "LBXEOPCT"],
    # 维生素
    "LBXVID": ["LBXVID", "LBXVIDMS"],
    "LBDVID": ["LBDVID", "LBDVIDLC"],
    # 肝功能 SI列名
    "LBXSALTI": ["LBXSALTI", "LBXSALTIS", "LBXSASSI"],
    "LBXSAT": ["LBXSAT", "LBXSATSI"],
    "LBXSGT": ["LBXSGT", "LBXSGTSI"],
    "LBXSAP": ["LBXSAP", "LBXSAPSI"],
    "LBXSTB": ["LBXSTB", "LBDSTBSI"],
    "LBXSAL": ["LBXSAL", "LBDSALSI"],
    # MCQ 疾病史 - 跨周期别名
    "MCQ010": ["MCQ010", "MCQ010A"],
    "MCQ035": ["MCQ035", "MCQ035A"],
    "MCQ053": ["MCQ053", "MCQ050"],
    "MCQ080": ["MCQ080", "MCQ080A"],
    "MCQ092": ["MCQ092", "MCQ090"],
    "MCQ120": ["MCQ120", "MCQ121"],
    "KIQ022": ["KIQ022", "KIQ020"],
    "KIQ025": ["KIQ025", "KIQ023"],
    "OSQ060": ["OSQ060", "OSQ050"],
}

# ============================================================================
# 预定义变量组
# ============================================================================
# 所有标签格式：中文名(英文缩写, 单位)
# 排序：核心人口学→体格→血压→临床检验→生活方式→详细人口学
VARIABLE_GROUPS = [
    # ===== 人口学(核心) =====
    {
        "group_name": "� 人口学(核心)",
        "group_key": "demo_core",
        "table_name": "DEMO",
        "description": "性别、年龄、种族、教育、收入、婚姻",
        "variables": [
            {"var_name": "SEQN", "var_label": "序号(SEQN)", "is_id": True},
            {"var_name": "RIAGENDR", "var_label": "性别(Gender)", "mapping": {1: "男", 2: "女"}},
            {"var_name": "RIDAGEYR", "var_label": "年龄(Age,岁)"},
            {"var_name": "RIDRETH1", "var_label": "种族(Race)",
             "mapping": {1: "墨西哥裔", 2: "其他西班牙裔", 3: "非西班牙裔白人", 4: "非西班牙裔黑人", 5: "其他种族"}},
            {"var_name": "DMDEDUC2", "var_label": "教育程度(Education)",
             "mapping": {1: "低于9年级", 2: "9-11年级", 3: "高中/同等学历", 4: "大学部分学分", 5: "大学及以上"}},
            {"var_name": "INDFMPIR", "var_label": "收入贫困比(PIR)"},
            {"var_name": "DMDMARTL", "var_label": "婚姻状况(Marital)",
             "mapping": {1: "已婚", 2: "丧偶", 3: "离婚", 4: "分居", 5: "从未结婚", 6: "同居"}},
            {"var_name": "SDMVPSU", "var_label": "PSU(SDMVPSU)"},
            {"var_name": "SDMVSTRA", "var_label": "分层(SDMVSTRA)"},
            {"var_name": "WTMEC2YR", "var_label": "体检权重(WT-MEC)"},
            {"var_name": "WTINT2YR", "var_label": "访谈权重(WT-Int)"},
        ],
    },
    # ===== 体格测量 =====
    {
        "group_name": "🩺 体格测量",
        "group_key": "bmx",
        "table_name": "BMX",
        "description": "身高、体重、BMI、腰围",
        "variables": [
            {"var_name": "BMXWT", "var_label": "体重(WT,kg)"},
            {"var_name": "BMXHT", "var_label": "身高(HT,cm)"},
            {"var_name": "BMXBMI", "var_label": "BMI(kg/m²)"},
            {"var_name": "BMXWAIST", "var_label": "腰围(WC,cm)"},
            {"var_name": "BMXHIP", "var_label": "臀围(HC,cm)"},
            {"var_name": "BMXLEG", "var_label": "腿长(Leg,cm)"},
            {"var_name": "BMXARML", "var_label": "臂长(Arm,cm)"},
            {"var_name": "BMXARMC", "var_label": "上臂围(AC,cm)"},
            {"var_name": "BMXTRI", "var_label": "三头肌皮褶(Triceps,mm)"},
            {"var_name": "BMXHEAD", "var_label": "头围(Head,cm)"},
        ],
    },
    # ===== 血压 =====
    {
        "group_name": "❤️ 血压",
        "group_key": "bpx",
        "table_name": "BPX",
        "description": "收缩压/舒张压、脉率",
        "variables": [
            {"var_name": "BPXSY1", "var_label": "收缩压-1(SBP1,mmHg)"},
            {"var_name": "BPXDI1", "var_label": "舒张压-1(DBP1,mmHg)"},
            {"var_name": "BPXSY2", "var_label": "收缩压-2(SBP2,mmHg)"},
            {"var_name": "BPXDI2", "var_label": "舒张压-2(DBP2,mmHg)"},
            {"var_name": "BPXSY3", "var_label": "收缩压-3(SBP3,mmHg)"},
            {"var_name": "BPXDI3", "var_label": "舒张压-3(DBP3,mmHg)"},
            {"var_name": "BPXSY4", "var_label": "收缩压-4(SBP4,mmHg)"},
            {"var_name": "BPXDI4", "var_label": "舒张压-4(DBP4,mmHg)"},
            {"var_name": "BPXPLS", "var_label": "脉率(PR,bpm)"},
        ],
    },
    # ===== 血压病史 =====
    {
        "group_name": "💊 高血压病史",
        "group_key": "bpq",
        "table_name": "BPQ",
        "description": "高血压诊断、用药",
        "variables": [
            {"var_name": "BPQ020", "var_label": "曾诊断高血压(HTN-Dx)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "BPQ030", "var_label": "服用降压药(HTN-Meds)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "BPQ080", "var_label": "医生告知高血脂(HLD-Dx)", "mapping": {1: "是", 2: "否"}},
        ],
    },
    # ===== 血脂 =====
    {
        "group_name": "🧪 血脂",
        "group_key": "lipid",
        "table_name": "TCHOL",
        "description": "总胆固醇、HDL、LDL、甘油三酯",
        "variables": [
            {"var_name": "LBXTC", "var_label": "总胆固醇(TC,mg/dL)"},
            {"var_name": "LBDHDL", "var_label": "高密度脂蛋白(HDL,mg/dL)"},
            {"var_name": "LBDLDL", "var_label": "低密度脂蛋白(LDL,mg/dL)"},
            {"var_name": "LBXTR", "var_label": "甘油三酯(TG,mg/dL)"},
        ],
    },
    # ===== 血糖/胰岛素 =====
    {
        "group_name": "🧪 血糖/胰岛素",
        "group_key": "glu",
        "table_name": "GLU",
        "description": "血糖、糖化血红蛋白、胰岛素、C肽",
        "variables": [
            {"var_name": "LBXGLU", "var_label": "空腹血糖(Glu,mg/dL)"},
            {"var_name": "LBXGH", "var_label": "糖化血红蛋白(HbA1c,%)"},
            {"var_name": "LBXIN", "var_label": "胰岛素(Ins,uU/mL)"},
            {"var_name": "LBXCP", "var_label": "C肽(C-Pep,nmol/L)"},
        ],
    },
    # ===== 甲状腺功能 =====
    {
        "group_name": "🧪 甲状腺功能",
        "group_key": "thy",
        "table_name": "THYROD",
        "description": "⚠️ 仅2007-2012年有数据 TSH/T4/T3/抗体",
        "variables": [
            {"var_name": "LBXTSH1", "var_label": "促甲状腺激素(TSH,mIU/L)"},
            {"var_name": "LBXT4", "var_label": "总甲状腺素(TT4,μg/dL)"},
            {"var_name": "LBXT4F", "var_label": "游离甲状腺素(FT4,ng/dL)"},
            {"var_name": "LBXT3F", "var_label": "游离三碘甲腺原氨酸(FT3,pg/mL)"},
            {"var_name": "LBXTT3", "var_label": "总三碘甲腺原氨酸(TT3,ng/dL)"},
            {"var_name": "LBXTGN", "var_label": "甲状腺球蛋白(Tg,ng/mL)"},
            {"var_name": "LBXTPO", "var_label": "TPO抗体(TPOAb,IU/mL)"},
            {"var_name": "LBXATG", "var_label": "Tg抗体(TgAb,IU/mL)"},
            {"var_name": "WTSA2YR", "var_label": "甲状腺权重(WT-Subset)"},
        ],
    },
    # ===== 肾功能 =====
    {
        "group_name": "🧪 肾功能",
        "group_key": "kidney",
        "table_name": "BIOPRO",
        "description": "肌酐、尿素氮、尿酸",
        "variables": [
            {"var_name": "LBXSCR", "var_label": "肌酐(Cr,mg/dL)"},
            {"var_name": "LBXSBU", "var_label": "尿素氮(BUN,mg/dL)"},
            {"var_name": "LBXSUA", "var_label": "尿酸(UA,mg/dL)"},
        ],
    },
    # ===== 肝功能 =====
    {
        "group_name": "🧪 肝功能",
        "group_key": "liver",
        "table_name": "BIOPRO",
        "description": "ALT、AST、GGT、胆红素、白蛋白",
        "variables": [
            {"var_name": "LBXSALTI", "var_label": "谷丙转氨酶(ALT,IU/L)"},
            {"var_name": "LBXSAT", "var_label": "谷草转氨酶(AST,IU/L)"},
            {"var_name": "LBXSGT", "var_label": "γ-谷氨酰转移酶(GGT,IU/L)"},
            {"var_name": "LBXSTB", "var_label": "总胆红素(TBIL,mg/dL)"},
            {"var_name": "LBXSAP", "var_label": "碱性磷酸酶(ALP,IU/L)"},
            {"var_name": "LBXSAL", "var_label": "白蛋白(ALB,g/dL)"},
        ],
    },
    # ===== 维生素 =====
    {
        "group_name": "🧪 维生素",
        "group_key": "vit",
        "table_name": "VID",
        "description": "维生素D、B12、叶酸",
        "variables": [
            {"var_name": "LBXVID", "var_label": "维生素D(VitD,nmol/L)"},
            {"var_name": "LBDVID", "var_label": "维生素D等级(VitD-Cat)",
             "mapping": {1: "充足", 2: "不足", 3: "缺乏"}},
            {"var_name": "LBXB12", "var_label": "维生素B12(VitB12,pg/mL)"},
            {"var_name": "LBXFOL", "var_label": "叶酸(Folate,ng/mL)"},
        ],
    },
    # ===== 炎症指标 =====
    {
        "group_name": "🧪 炎症/血常规",
        "group_key": "inf",
        "table_name": "CBC",
        "description": "白细胞、红细胞、血小板、CRP",
        "variables": [
            {"var_name": "LBXWBC", "var_label": "白细胞(WBC,×10⁹/L)"},
            {"var_name": "LBXNEUT", "var_label": "中性粒细胞(NEUT,×10⁹/L)"},
            {"var_name": "LBXLY", "var_label": "淋巴细胞(LY,×10⁹/L)"},
            {"var_name": "LBXMONO", "var_label": "单核细胞(MONO,×10⁹/L)"},
            {"var_name": "LBXEOS", "var_label": "嗜酸粒细胞(EO,×10⁹/L)"},
            {"var_name": "LBXRBC", "var_label": "红细胞(RBC,×10¹²/L)"},
            {"var_name": "LBXHGB", "var_label": "血红蛋白(Hb,g/dL)"},
            {"var_name": "LBXHCT", "var_label": "红细胞压积(HCT,%)"},
            {"var_name": "LBXPLT", "var_label": "血小板(PLT,×10⁹/L)"},
            {"var_name": "LBXCRP", "var_label": "超敏C反应蛋白(hs-CRP,mg/dL)"},
        ],
    },
    # ===== 重金属 =====
    {
        "group_name": "🧪 重金属",
        "group_key": "metal",
        "table_name": "PBCD",
        "description": "血铅、镉、汞、硒",
        "variables": [
            {"var_name": "LBXBPB", "var_label": "血铅(Pb,μg/dL)"},
            {"var_name": "LBXBCD", "var_label": "血镉(Cd,μg/L)"},
            {"var_name": "LBXTHG", "var_label": "总汞(Hg,μg/L)"},
            {"var_name": "LBXBSE", "var_label": "血硒(Se,μg/L)"},
        ],
    },
    # ===== 性激素 =====
    {
        "group_name": "🧪 性激素",
        "group_key": "sexhrm",
        "table_name": "SEXHRM",
        "description": "睾酮、雌二醇、SHBG",
        "variables": [
            {"var_name": "LBXTST", "var_label": "总睾酮(T,ng/dL)"},
            {"var_name": "LBXEST", "var_label": "雌二醇(E2,pg/mL)"},
            {"var_name": "LBXSHBG", "var_label": "性激素结合球蛋白(SHBG,nmol/L)"},
        ],
    },
    # ===== 吸烟 =====
    {
        "group_name": "🚬 吸烟",
        "group_key": "smq",
        "table_name": "SMQ",
        "description": "吸烟史、目前吸烟量",
        "variables": [
            {"var_name": "SMQ020", "var_label": "曾吸100支烟(Smok100+)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "SMQ040", "var_label": "目前吸烟频率(SmokNow)",
             "mapping": {1: "每天吸", 2: "偶尔吸", 3: "不吸"}},
            {"var_name": "SMQ050Q", "var_label": "日均吸烟量(Cig/d,支)"},
        ],
    },
    # ===== 饮酒 =====
    {
        "group_name": "🍺 饮酒",
        "group_key": "alq",
        "table_name": "ALQ",
        "description": "饮酒史、频率、暴饮",
        "variables": [
            {"var_name": "ALQ101", "var_label": "曾饮酒12次以上(EverDrink)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "ALQ110", "var_label": "过去1年饮酒(Drink1y)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "ALQ120Q", "var_label": "饮酒频率(Freq,d/月)"},
            {"var_name": "ALQ130", "var_label": "日均饮酒量(Drinks/d,杯)"},
        ],
    },
    # ===== 糖尿病 =====
    {
        "group_name": "💉 糖尿病",
        "group_key": "diq",
        "table_name": "DIQ",
        "description": "糖尿病诊断、治疗",
        "variables": [
            {"var_name": "DIQ010", "var_label": "糖尿病诊断(DM-Dx)",
             "mapping": {1: "是", 2: "否", 3: "临界/妊娠期糖尿病"}},
            {"var_name": "DIQ050", "var_label": "使用胰岛素(Insulin)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "DIQ070", "var_label": "口服降糖药(OHA)", "mapping": {1: "是", 2: "否"}},
        ],
    },
    # ===== 体力活动 =====
    {
        "group_name": "🏃 体力活动",
        "group_key": "paq",
        "table_name": "PAQ",
        "description": "工作/休闲活动、久坐",
        "variables": [
            {"var_name": "PAQ605", "var_label": "剧烈工作活动(VigWork)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "PAQ620", "var_label": "中等工作活动(ModWork)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "PAQ650", "var_label": "休闲步行/骑车(WalkBike)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "PAQ665", "var_label": "每日静坐(Sedentary,h)"},
        ],
    },
    # ===== 睡眠 =====
    {
        "group_name": "😴 睡眠",
        "group_key": "slq",
        "table_name": "SLQ",
        "description": "睡眠时长、障碍、打鼾",
        "variables": [
            {"var_name": "SLD010H", "var_label": "睡眠时长(Sleep,h/d)"},
            {"var_name": "SLQ050", "var_label": "睡眠障碍(Sleep-Dx)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "SLQ120", "var_label": "打鼾(Snore)",
             "mapping": {0: "从不", 1: "很少", 2: "有时", 3: "经常"}},
        ],
    },
    # ===== 心理健康 =====
    {
        "group_name": "😊 心理健康(PHQ-9)",
        "group_key": "dpq",
        "table_name": "DPQ",
        "description": "PHQ-9抑郁症状",
        "variables": [
            {"var_name": "DPQ010", "var_label": "兴趣丧失(Anhedonia)",
             "mapping": {0: "无", 1: "几天", 2: ">一半", 3: "几乎每天"}},
            {"var_name": "DPQ020", "var_label": "情绪低落(Depressed)",
             "mapping": {0: "无", 1: "几天", 2: ">一半", 3: "几乎每天"}},
            {"var_name": "DPQ030", "var_label": "睡眠障碍(Sleep-Prob)",
             "mapping": {0: "无", 1: "几天", 2: ">一半", 3: "几乎每天"}},
            {"var_name": "DPQ040", "var_label": "疲乏无力(Fatigue)",
             "mapping": {0: "无", 1: "几天", 2: ">一半", 3: "几乎每天"}},
            {"var_name": "DPQ090", "var_label": "自伤念头(Suicidal)",
             "mapping": {0: "无", 1: "几天", 2: ">一半", 3: "几乎每天"}},
        ],
    },
    # ===== 医疗可及性 =====
    {
        "group_name": "🏥 医疗可及性",
        "group_key": "huq",
        "table_name": "HUQ",
        "description": "自评健康、医保、就医",
        "variables": [
            {"var_name": "HUQ010", "var_label": "自评健康(Health-Gen)",
             "mapping": {1: "非常好", 2: "很好", 3: "好", 4: "一般", 5: "差"}},
            {"var_name": "HUQ050", "var_label": "有医疗保险(Insurance)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "HUQ030", "var_label": "有固定就医点(RegularCare)", "mapping": {1: "是", 2: "否"}},
        ],
    },
    # ===== 处方药 =====
    {
        "group_name": "💊 处方药",
        "group_key": "rxq",
        "table_name": "RXQ_RX",
        "description": "处方药使用数量",
        "variables": [
            {"var_name": "RXDUSE", "var_label": "使用处方药(Rx-Use)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "RXDCOUNT", "var_label": "处方药数(Rx-Count,n)"},
        ],
    },
    # ===== 既往病史(MCQ) =====
    {
        "group_name": "🏥 既往病史(自报疾病)",
        "group_key": "mcq",
        "table_name": "MCQ",
        "description": "医生曾告知的疾病：心脏病、肝病、肾病、关节炎、癌症等",
        "variables": [
            {"var_name": "MCQ010", "var_label": "哮喘(Asthma)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ035", "var_label": "关节炎(Arthritis)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ053", "var_label": "心力衰竭(CHF)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ055", "var_label": "冠心病(CHD)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ057", "var_label": "心绞痛(Angina)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ059", "var_label": "心肌梗死(HeartAttack)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ062", "var_label": "中风(Stroke)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ065", "var_label": "肺气肿(Emphysema)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ067", "var_label": "慢性支气管炎(ChrBronch)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ070", "var_label": "COPD", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ080", "var_label": "癌症(Cancer)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ084", "var_label": "痴呆/阿尔茨海默(Dementia)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ092", "var_label": "甲状腺问题(Thyroid)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ120", "var_label": "慢性肝病(Liver)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ160A", "var_label": "骨关节炎(Osteoarth)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ160B", "var_label": "类风湿关节炎(RA)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ160F", "var_label": "痛风(Gout)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "MCQ170M", "var_label": "甲状腺病(ThyroidDz)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "KIQ022", "var_label": "肾衰竭(KidneyFail)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "KIQ025", "var_label": "慢性肾病(CKD)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "OSQ060", "var_label": "骨质疏松(Osteoporosis)", "mapping": {1: "是", 2: "否"}},
        ],
    },
    # ===== 饮食 =====
    {
        "group_name": "🍽️ 饮食",
        "group_key": "diet",
        "table_name": "DR1TOT",
        "description": "24h膳食回顾(总能量、宏量营养素)",
        "variables": [
            {"var_name": "DR1TKCAL", "var_label": "总能量(Energy,kcal)"},
            {"var_name": "DR1TPROT", "var_label": "蛋白质(Protein,g)"},
            {"var_name": "DR1TCARB", "var_label": "碳水化合物(CHOs,g)"},
            {"var_name": "DR1TTFAT", "var_label": "总脂肪(TFAT,g)"},
            {"var_name": "DR1TSFAT", "var_label": "饱和脂肪(SFAT,g)"},
            {"var_name": "DR1TFIBE", "var_label": "膳食纤维(Fiber,g)"},
            {"var_name": "DR1TSODI", "var_label": "钠(Na,mg)"},
            {"var_name": "DR1TALCO", "var_label": "酒精(Alcohol,g)"},
        ],
    },
    # ===== 人口学(详细) =====
    {
        "group_name": "📋 人口学(详细)",
        "group_key": "demo_detail",
        "table_name": "DEMO",
        "description": "年龄月龄、种族亚裔、教育青少年、收入、出生地等",
        "variables": [
            {"var_name": "RIDAGEMN", "var_label": "年龄(Age,月)"},
            {"var_name": "RIDRETH3", "var_label": "种族-含亚裔(Race-Asian)",
             "mapping": {1: "墨西哥裔", 2: "其他西班牙裔", 3: "非西班牙裔白人",
                         4: "非西班牙裔黑人", 6: "非西班牙裔亚裔", 7: "其他"}},
            {"var_name": "DMDEDUC3", "var_label": "教育-青少年版(Edu-Youth)",
             "mapping": {0: "未上学", 1: "1-5年级", 2: "6-8年级", 3: "9-11年级", 4: "12年级", 5: "大学以上"}},
            {"var_name": "INDHHIN2", "var_label": "家庭年收入(HHIncome,USD)"},
            {"var_name": "DMDBORN4", "var_label": "出生国(Birth-US)", "mapping": {1: "美国", 2: "其他"}},
            {"var_name": "DMDCITZN", "var_label": "公民(Citizen)", "mapping": {1: "是", 2: "否"}},
            {"var_name": "WTINT2YR", "var_label": "访谈权重(WT-Int)"},
            {"var_name": "WTMEC2YR", "var_label": "体检权重(WT-MEC)"},
            {"var_name": "WTSA2YR", "var_label": "甲状腺权重(WT-Subset)"},
            {"var_name": "SDDSRVYR", "var_label": "调查周期编号(Cycle-N)"},
            {"var_name": "RIDSTATR", "var_label": "参与状态(Status)",
             "mapping": {1: "访谈+体检都完成", 2: "仅完成访谈"}},
            {"var_name": "RIDEXPRG", "var_label": "怀孕状态(Pregnant)",
             "mapping": {1: "是", 2: "否", 3: "不确定"}},
        ],
    },
    # ===== 口腔健康 =====
    {
        "group_name": "🦷 口腔健康",
        "group_key": "ohx",
        "table_name": "OHX",
        "description": "龋齿、补牙、缺失",
        "variables": [
            {"var_name": "OHX02HC", "var_label": "龋齿数(Caries,n)"},
            {"var_name": "OHX04HC", "var_label": "补牙数(Filling,n)"},
            {"var_name": "OHX06HC", "var_label": "缺失数(Missing,n)"},
        ],
    },
]

# ============================================================================
# 死亡率数据定义
# ============================================================================
# CDC Public-Use Linked Mortality Files 下载地址：
# https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/
#
# 固定宽度格式解析说明:
#   SEQN: 1-6
#   ELIGSTAT: 15
#   MORTSTAT: 16
#   UCOD_LEADING: 17-19
#   DIABETES: 20
#   HYPERTEN: 21
#   PERMTH_INT: 43-45
#   PERMTH_EXM: 46-48

MORTALITY_FWF_SPEC = [
    ("SEQN", 1, 6, "i"),
    ("ELIGSTAT", 15, 15, "i"),
    ("MORTSTAT", 16, 16, "i"),
    ("UCOD_LEADING", 17, 19, "i"),
    ("DIABETES", 20, 20, "i"),
    ("HYPERTEN", 21, 21, "i"),
    ("PERMTH_INT", 43, 45, "i"),
    ("PERMTH_EXM", 46, 48, "i"),
]

MORTALITY_GROUPS = [
    {
        "group_name": "💀 死亡结局(CDC链接数据)",
        "group_key": "mort",
        "table_name": "MORT",
        "description": "死亡率链接数据：生存状态、随访时间、死因",
        "variables": [
            {"var_name": "ELIGSTAT", "var_label": "死亡率随访资格",
             "mapping": {1: "合格", 2: "<18岁(不公开)", 3: "不合格"}},
            {"var_name": "MORTSTAT", "var_label": "最终死亡状态",
             "mapping": {0: "存活", 1: "死亡"}},
            {"var_name": "UCOD_LEADING", "var_label": "主要死因分类",
             "mapping": {1: "心脏病", 2: "恶性肿瘤", 3: "慢性下呼吸道疾病",
                         4: "意外伤害", 5: "脑血管疾病", 6: "阿尔茨海默病",
                         7: "糖尿病", 8: "流感/肺炎", 9: "肾炎/肾病", 10: "其他"}},
            {"var_name": "DIABETES", "var_label": "死因-DM标志",
             "mapping": {0: "否", 1: "是"}},
            {"var_name": "HYPERTEN", "var_label": "死因-HTN标志",
             "mapping": {0: "否", 1: "是"}},
            {"var_name": "PERMTH_INT", "var_label": "随访时间(月,从访谈算起)"},
            {"var_name": "PERMTH_EXM", "var_label": "随访时间(月,从体检算起)"},
        ],
    },
]

# ============================================================================
# 标签映射
# ============================================================================

def build_var_label_map():
    m = {"SEQN": "序号(SEQN)"}
    for g in VARIABLE_GROUPS:
        for v in g["variables"]:
            m[v["var_name"]] = v["var_label"]
    for g in MORTALITY_GROUPS:
        for v in g["variables"]:
            m[v["var_name"]] = v["var_label"]
    return m

VAR_LABEL_MAP = build_var_label_map()

def get_var_label(n):
    return VAR_LABEL_MAP.get(n, n)


# ============================================================================
# 下载器
# ============================================================================

class Downloader:
    def __init__(self):
        self.cache_dir = CACHE_DIR
        os.makedirs(self.cache_dir, exist_ok=True)
        self.ssl_context = ssl.create_default_context()
        self.ssl_context.check_hostname = False
        self.ssl_context.verify_mode = ssl.CERT_NONE

    def get_xpt_url(self, table_base, cycle_suffix, start_year):
        if cycle_suffix == "A":
            xpt_filename = f"{table_base}.xpt"
        else:
            xpt_filename = f"{table_base}_{cycle_suffix}.xpt"
        url = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{start_year}/DataFiles/{xpt_filename}"
        return url, xpt_filename

    def download_xpt(self, table_base, cycle_suffix, start_year, callback=None):
        url, fname = self.get_xpt_url(table_base, cycle_suffix, start_year)
        cache_path = os.path.join(self.cache_dir, fname)

        if os.path.exists(cache_path) and os.path.getsize(cache_path) > 1024:
            if callback: callback(f"缓存命中: {fname}")
            return cache_path, True

        try:
            if callback: callback(f"下载中: {url}")
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, context=self.ssl_context, timeout=120) as resp:
                data = resp.read()
            if data.strip().startswith(b"<"):
                if callback: callback(f"文件不存在: {fname}")
                return None, False
            with open(cache_path, "wb") as f:
                f.write(data)
            if callback: callback(f"下载完成: {fname} ({len(data)//1024}KB)")
            return cache_path, True
        except urllib.error.HTTPError as e:
            if callback: callback(f"HTTP {e.code}: {fname}")
            return None, False
        except Exception as e:
            if callback: callback(f"下载失败 {fname}: {e}")
            return None, False

    def download_mortality(self, cycle_suffix, start_year, callback=None):
        """下载死亡率数据(.dat固定宽度文件)"""
        cdc_year_map = {
            "A": ("1999_2000", "NHANES_1999_2000_MORT_2019_PUBLIC.dat"),
            "B": ("2001_2002", "NHANES_2001_2002_MORT_2019_PUBLIC.dat"),
            "C": ("2003_2004", "NHANES_2003_2004_MORT_2019_PUBLIC.dat"),
            "D": ("2005_2006", "NHANES_2005_2006_MORT_2019_PUBLIC.dat"),
            "E": ("2007_2008", "NHANES_2007_2008_MORT_2019_PUBLIC.dat"),
            "F": ("2009_2010", "NHANES_2009_2010_MORT_2019_PUBLIC.dat"),
            "G": ("2011_2012", "NHANES_2011_2012_MORT_2019_PUBLIC.dat"),
            "H": ("2013_2014", "NHANES_2013_2014_MORT_2019_PUBLIC.dat"),
            "I": ("2015_2016", "NHANES_2015_2016_MORT_2019_PUBLIC.dat"),
            "J": ("2017_2018", "NHANES_2017_2018_MORT_2019_PUBLIC.dat"),
        }
        if cycle_suffix not in cdc_year_map:
            if callback: callback(f"此周期无死亡率数据: {cycle_suffix}")
            return None, False

        year_tag, dat_filename = cdc_year_map[cycle_suffix]
        url = f"https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/{dat_filename}"
        cache_path = os.path.join(self.cache_dir, dat_filename)

        if os.path.exists(cache_path) and os.path.getsize(cache_path) > 100:
            if callback: callback(f"死亡率数据已缓存: {dat_filename}")
            return cache_path, True

        try:
            if callback: callback(f"下载死亡率数据: {url}")
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, context=self.ssl_context, timeout=120) as resp:
                data = resp.read()
            with open(cache_path, "wb") as f:
                f.write(data)
            if callback: callback(f"死亡率数据下载完成: {dat_filename} ({len(data)//1024}KB)")
            return cache_path, True
        except Exception as e:
            if callback: callback(f"死亡率数据下载失败: {e}")
            return None, False


# ============================================================================
# XPT 读取
# ============================================================================

def read_xpt(file_path):
    try:
        return pd.read_sas(file_path, format="xport", encoding="utf-8")
    except:
        try:
            return pd.read_sas(file_path, format="xport", encoding="latin-1")
        except Exception as e:
            print(f"XPT读取失败 {file_path}: {e}")
            return None


# ============================================================================
# 死亡率数据解析(固定宽度)
# ============================================================================

def read_mortality(file_path):
    """解析死亡率 .dat 固定宽度文件（逐行读取，支持CRLF/LF换行）"""
    if not file_path or not os.path.exists(file_path):
        return None
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            data = []
            with open(file_path, "r", encoding="ascii", errors="replace") as f:
                lines = f.readlines()

            for line in lines:
                raw_line = line.rstrip('\n\r')
                if len(raw_line) < 21:  # 至少需要读到HYPERTEN(位置21)
                    continue
                try:
                    seqn_str = raw_line[0:6].strip()
                    seqn = int(seqn_str) if seqn_str.isdigit() else -1
                except:
                    seqn = -1
                if seqn <= 0:
                    continue

                def get_int(start, end):
                    """从1-indexed位置start到end取值"""
                    try:
                        val = raw_line[start-1:end].strip()
                        return int(val) if val else -1
                    except:
                        return -1

                data.append({
                    "SEQN": seqn,
                    "ELIGSTAT": get_int(15, 15),
                    "MORTSTAT": get_int(16, 16),
                    "UCOD_LEADING": get_int(17, 19),
                    "DIABETES": get_int(20, 20),
                    "HYPERTEN": get_int(21, 21),
                    "PERMTH_INT": get_int(43, 45),
                    "PERMTH_EXM": get_int(46, 48),
                })

        if not data:
            return None
        df = pd.DataFrame(data)
        # SEQN保持为Int64，其他字段用float64（NaN兼容）
        df["SEQN"] = pd.to_numeric(df["SEQN"], errors='coerce').astype('Int64')
        for c in df.columns:
            if c == "SEQN":
                continue
            df[c] = df[c].replace(-1, np.nan).astype("float64")
        
        # CDC规定：ELIGSTAT=1为死亡率合格样本，=2为年龄<18岁不公开，=3为其他不合格
        # 只保留合格样本，确保样本量与CDC官方一致
        before = len(df)
        df = df[df["ELIGSTAT"] == 1].copy()
        if len(df) < before:
            pass  # 过滤掉 {before - len(df)} 条不合格记录
        
        return df
    except Exception as e:
        print(f"死亡率文件解析失败: {e}")
        return None


# ============================================================================
# 数据清理
# ============================================================================

MISSING_VALUES = [
    77777, 77777.0, 99999, 99999.0, 88888, 88888.0,
    7777, 7777.0, 9999, 9999.0, 8888, 8888.0,
    777, 777.0, 999, 999.0, 888, 888.0,
    77, 77.0, 99, 99.0, 88, 88.0,
    # 注意: 7,8,9 不再作为通用缺失值，因为它们是睡眠时长等变量的正常值
]

# 睡眠时长等变量的正常值范围（这些值不应被视为缺失值）
NORMAL_VALUE_RANGES = {
    "sleep": {"var_names": ["SLD010H", "睡眠时长"], "min": 1, "max": 24},
    "snore": {"var_names": ["SLQ120", "打鼾"], "min": 0, "max": 3},
}

def apply_missing_clean(series, var_label=None, var_name=None):
    """将NHANES缺失值标记替换为NaN
    
    Args:
        series: 数据系列
        var_label: 变量标签（中文），用于判断是否为特殊变量
        var_name: 变量名（英文），用于判断是否为特殊变量
    """
    s = series.copy()
    if s.dtype not in [np.float64, np.float32, np.int64, np.int32]:
        return s
    
    # 检查是否为睡眠时长等特殊变量
    # 优先使用var_name（原始XPT列名，更可靠），避免var_label中的中文关键词被误匹配
    is_sleep_var = False
    if var_name:
        for key, cfg in NORMAL_VALUE_RANGES.items():
            if any(vn.upper() == (var_name or "").upper() for vn in cfg["var_names"]):
                is_sleep_var = True
                min_val = cfg["min"]
                max_val = cfg["max"]
                s = s.apply(lambda x: np.nan if pd.isna(x) or x < min_val or x > max_val else x)
                return s
    # var_label作为兜底（仅在var_name为空时使用）
    if not is_sleep_var and var_label:
        for key, cfg in NORMAL_VALUE_RANGES.items():
            if any(vn in var_label for vn in cfg["var_names"]):
                min_val = cfg["min"]
                max_val = cfg["max"]
                s = s.apply(lambda x: np.nan if pd.isna(x) or x < min_val or x > max_val else x)
                return s
    
    # 非特殊变量，使用通用缺失值替换
    s = s.replace(MISSING_VALUES, np.nan)
    
    # SAS特殊缺失值处理：pandas读取XPT时，SAS的缺失值编码（., .A, .B等）
    # 会变成极小的浮点数（如 5.397605e-79 = 2^(-260)），
    # 这些值不在MISSING_VALUES里但也不是有效数据。
    # 特征：绝对值在 1e-20 ~ 1e-100 之间且不是0
    if s.dtype in [np.float64, np.float32]:
        tiny_mask = (s.abs() > 1e-100) & (s.abs() < 1e-15)
        if tiny_mask.any():
            n_tiny = tiny_mask.sum()
            if n_tiny > 0:
                s = s.mask(tiny_mask, np.nan)
    
    return s


# ============================================================================
# 单位换算（美制 → 中国常用单位）
# ============================================================================
# 格式: (列名片段, 换算系数, 目标单位, 精确度)
# 规则：若列名包含关键词则应用换算
UNIT_CONVERSIONS = [
    # 血脂
    ("总胆固醇(TC,mg/dL)", 0.02586, "总胆固醇(TC,mmol/L)", 2),
    ("高密度脂蛋白(HDL,mg/dL)", 0.02586, "高密度脂蛋白(HDL,mmol/L)", 2),
    ("低密度脂蛋白(LDL,mg/dL)", 0.02586, "低密度脂蛋白(LDL,mmol/L)", 2),
    ("甘油三酯(TG,mg/dL)", 0.01129, "甘油三酯(TG,mmol/L)", 2),
    # 血糖
    ("空腹血糖(Glu,mg/dL)", 0.0555, "空腹血糖(Glu,mmol/L)", 1),
    # 肾功能（国内用μmol/L）
    ("肌酐(Cr,mg/dL)", 88.4, "肌酐(Cr,μmol/L)", 0),
    ("尿酸(UA,mg/dL)", 59.48, "尿酸(UA,μmol/L)", 0),
    ("尿素氮(BUN,mg/dL)", 0.357, "尿素氮(BUN,mmol/L)", 1),
    # 甲状腺
    ("促甲状腺激素(TSH,mIU/L)", 1.0, "促甲状腺激素(TSH,mIU/L)", 2),
    ("总甲状腺素(TT4,μg/dL)", 12.87, "总甲状腺素(TT4,nmol/L)", 1),
    ("游离甲状腺素(FT4,ng/dL)", 12.87, "游离甲状腺素(FT4,pmol/L)", 1),
    ("总三碘甲腺原氨酸(TT3,ng/dL)", 0.0154, "总三碘甲腺原氨酸(TT3,nmol/L)", 3),
    ("游离三碘甲腺原氨酸(FT3,pg/mL)", 1.536, "游离三碘甲腺原氨酸(FT3,pmol/L)", 1),
    # 肝功能
    ("总胆红素(TBIL,mg/dL)", 17.1, "总胆红素(TBIL,μmol/L)", 1),
    ("白蛋白(ALB,g/dL)", 10, "白蛋白(ALB,g/L)", 0),
    # 其他
    ("血铅(Pb,μg/dL)", 0.0483, "血铅(Pb,μmol/L)", 2),
    ("血红蛋白(Hb,g/dL)", 10, "血红蛋白(Hb,g/L)", 0),
]

def apply_unit_conversion(df, enabled=True):
    """对DataFrame应用单位换算，新增换算后的列"""
    if not enabled:
        return df
    df = df.copy()
    for keyword, factor, new_label, decimals in UNIT_CONVERSIONS:
        for col in df.columns:
            if keyword in col and col != new_label:
                try:
                    converted = pd.to_numeric(df[col], errors='coerce') * factor
                    df[new_label] = converted.round(decimals)
                except:
                    pass
    return df


COLUMN_ORDER_PREFIX = [
    "序号(SEQN)", "调查周期",
    "性别(Gender)", "年龄(Age,岁)", "种族(Race)", "教育程度(Education)", "收入贫困比(PIR)", "婚姻状况(Marital)",
    "体重(WT,kg)", "身高(HT,cm)", "BMI(kg/m²)", "腰围(WC,cm)",
    "收缩压-1(SBP1,mmHg)", "舒张压-1(DBP1,mmHg)",
    "促甲状腺激素(TSH,mIU/L)",
    "游离三碘甲腺原氨酸(FT3,pmol/L)", "游离甲状腺素(FT4,pmol/L)",
    "总胆固醇(TC,mmol/L)", "甘油三酯(TG,mmol/L)", "高密度脂蛋白(HDL,mmol/L)", "低密度脂蛋白(LDL,mmol/L)",
    "空腹血糖(Glu,mmol/L)",
    "谷丙转氨酶(ALT,IU/L)", "谷草转氨酶(AST,IU/L)", "γ-谷氨酰转移酶(GGT,IU/L)",
    "肌酐(Cr,μmol/L)", "尿酸(UA,μmol/L)",
]

def clean_and_reorder_columns(df, keep_labels=None):
    """去除美标原始列、非核心列，按逻辑顺序排列
    keep_labels: set of column labels to keep (来自用户所选变量组)
    """
    # 需去除的美标关键词
    US_UNITS = ["(TC,mg/dL)", "(HDL,mg/dL)", "(LDL,mg/dL)", "(TG,mg/dL)",
                "(Glu,mg/dL)", "(Cr,mg/dL)", "(UA,mg/dL)", "(BUN,mg/dL)",
                "(TT4,μg/dL)", "(FT4,ng/dL)", "(TT3,ng/dL)", "(FT3,pg/mL)",
                "(TBIL,mg/dL)", "(ALB,g/dL)", "(Hb,g/dL)", "(Pb,μg/dL)"]

    # 只保留的核心列模式（其他列如吸烟/饮酒/饮食/病史/多余体格测量等全部丢弃）
    CORE_PATTERNS = [
        "序号(SEQN)", "调查周期",
        "性别", "年龄", "种族", "Race", "教育程度", "收入贫困比", "婚姻状况",
        "体重(WT,kg)", "身高(HT,cm)", "BMI", "腰围",
        "收缩压", "舒张压", "脉率", "高血压分级",
        "TSH", "FT3", "FT4", "TT3", "TT4",
        "TC,mmol", "TG,mmol", "HDL,mmol", "LDL,mmol",
        "Glu,mmol", "空腹血糖",
        "ALT", "AST", "GGT", "谷丙", "谷草",
        "Cr,μmol", "UA,μmol", "尿酸", "肌酐", "尿素氮", "BUN",
        "TBIL", "胆红素", "白蛋白", "ALB",
        "血红蛋白", "Hb,g", "白细胞", "WBC", "血小板",
        "体检权重", "访谈权重", "WT-MEC", "WT-Int", "WTSA", "WT-Subset",
        "SDMVPSU", "SDMVSTRA", "PSU", "Strata",
        "死亡率", "MORT", "死亡", "随访资格",
        "TSH/FT4", "FT3/FT4", "TyG", "动脉粥样硬化", "非HDL", "HDL-C达标", "血脂达标",
        "降脂药", "甲状腺药物", "处方药数", "处方药(Rx",
        "使用处方药", "Rx-Use", "Rx-Count",
        "调查周期编号", "Cycle-N",
        "异常", "异常值标记",
        "睡眠", "Sleep", "打鼾", "Snore",
    ]
    # 明确要丢弃的关键词（即使被CORE_PATTERNS匹配到）
    EXPLICIT_DROP = ["哮喘", "Asthma", "腿长", "臂长", "上臂围", "头围",
                     "Leg", "Arm", "AC", "Head", "年龄(Age,月)", "Cycle-N"]

    drop_cols = []
    for col in df.columns:
        # 检查美标
        for us_kw in US_UNITS:
            if us_kw in col:
                cn_kw = (us_kw.replace("(TC,mg/dL)","(TC,mmol/L)").replace("(HDL,mg/dL)","(HDL,mmol/L)")
                         .replace("(LDL,mg/dL)","(LDL,mmol/L)").replace("(TG,mg/dL)","(TG,mmol/L)")
                         .replace("(Glu,mg/dL)","(Glu,mmol/L)").replace("(Cr,mg/dL)","(Cr,μmol/L)")
                         .replace("(UA,mg/dL)","(UA,μmol/L)").replace("(BUN,mg/dL)","(BUN,mmol/L)")
                         .replace("(TT4,μg/dL)","(TT4,nmol/L)").replace("(FT4,ng/dL)","(FT4,pmol/L)")
                         .replace("(TT3,ng/dL)","(TT3,nmol/L)").replace("(FT3,pg/mL)","(FT3,pmol/L)")
                         .replace("(TBIL,mg/dL)","(TBIL,μmol/L)").replace("(ALB,g/dL)","(ALB,g/L)")
                         .replace("(Hb,g/dL)","(Hb,g/L)").replace("(Pb,μg/dL)","(Pb,μmol/L)"))
                if any(cn_kw in c for c in df.columns):
                    drop_cols.append(col)
                break
        else:
            # 不是美标列，判读是否是核心列
            # 如果传入了 keep_labels（用户所选变量组的标签），优先匹配
            if keep_labels:
                # 检查列名是否匹配任何所选变量组的标签
                col_clean = col.upper().replace(" ","").replace("-","").replace("(","").replace(")","")
                is_kept = False
                for kl in keep_labels:
                    kl_clean = kl.upper().replace(" ","").replace("-","").replace("(","").replace(")","")
                    if kl_clean in col_clean or col_clean in kl_clean:
                        is_kept = True
                        break
                if not is_kept:
                    drop_cols.append(col)
            else:
                col_up = col.upper().replace(" ","").replace("-","").replace("(","").replace(")","")
                is_core = any(p.upper().replace(" ","").replace("-","").replace("(","").replace(")","") in col_up for p in CORE_PATTERNS)
                if not is_core:
                    drop_cols.append(col)

    # 确保不删除序号和周期
    safe_keep = {"序号(SEQN)", "调查周期"}
    drop_cols = [c for c in drop_cols if c not in safe_keep]

    # 添加显式丢弃列（仅在静态白名单模式下生效，动态白名单模式下保留所有用户选择的变量）
    if not keep_labels:
        for col in df.columns:
            for kw in EXPLICIT_DROP:
                if kw.upper() in col.upper().replace(" ","").replace("-",""):
                    if col not in drop_cols and col not in safe_keep:
                        drop_cols.append(col)
                    break

    if drop_cols:
        # 记录被丢弃的前20个列名
        dropped_preview = ', '.join(sorted(drop_cols[:20]))
        if len(drop_cols) > 20:
            dropped_preview += f' … 共{len(drop_cols)}列'
        print(f"  [清洗] 已丢弃非核心列: {dropped_preview}")
        df.drop(columns=[c for c in drop_cols if c in df.columns], inplace=True, errors="ignore")

    # 重新排列
    all_cols = list(df.columns)
    ordered = [c for c in COLUMN_ORDER_PREFIX if c in all_cols]
    remaining = [c for c in all_cols if c not in ordered and c not in ("序号(SEQN)", "调查周期")]
    final_order = ordered + remaining
    final_order = [c for c in final_order if c in all_cols]
    df = df[final_order]
    return df


# ============================================================================
# 变量探索器 - 获取某个周期所有可用变量
# ============================================================================

class VariableExplorer:
    """从已下载的XPT文件中提取所有可用变量名"""

    def __init__(self, downloader):
        self.downloader = downloader
        self.cache = {}  # (table_name, suffix) -> [column_names]

    def get_variables(self, table_name, cycle_suffix, start_year):
        """获取某个表在某个周期的所有变量名"""
        key = (table_name, cycle_suffix)
        if key in self.cache:
            return self.cache[key]

        fp, ok = self.downloader.download_xpt(table_name, cycle_suffix, start_year)
        if not ok or not fp:
            return []

        df = read_xpt(fp)
        if df is None:
            return []

        cols = list(df.columns)
        self.cache[key] = cols
        return cols

    def get_all_tables_variables(self, cycle_suffix, start_year, table_names):
        """获取多个表的所有变量"""
        result = {}
        for tn in table_names:
            cols = self.get_variables(tn, cycle_suffix, start_year)
            if cols:
                result[tn] = cols
        return result


# ============================================================================
# 主处理引擎
# ============================================================================

# ============================================================================
# 依赖检查
# ============================================================================

def check_dependencies(selected_cycle_suffixes, selected_group_keys, profile):
    """运行前检查配置矛盾，返回警告列表"""
    warnings = []
    if profile is None:
        return warnings
    cfg = profile.config
    med_classes = cfg["exclusion"].get("medication_classes", [])
    if med_classes:
        if "rxq" not in selected_group_keys:
            warnings.append("⚠️ 需要药物排除但未勾选「处方药(RXQ_RX)」组，药物排除将失败！")
    disease_classes = cfg["exclusion"].get("disease_classes", [])
    if disease_classes:
        needed_mcq_vars = set()
        for dc in disease_classes:
            if dc in DISEASE_EXCLUSION_CATEGORIES:
                for v in DISEASE_EXCLUSION_CATEGORIES[dc]["vars"]:
                    needed_mcq_vars.add(v)
        if needed_mcq_vars and "mcq" not in selected_group_keys:
            warnings.append(f"⚠️ 需要排除疾病({', '.join(disease_classes)})但未勾选「既往病史(MCQ)」组，疾病排除将跳过！")
    return warnings


class NhanesEngine:
    def __init__(self):
        self.downloader = Downloader()
        self.explorer = VariableExplorer(self.downloader)
        self.log_callback = None
        self.is_running = False

    def log(self, msg):
        if self.log_callback:
            self.log_callback(msg)
        print(msg)

    def run(self, selected_cycle_suffixes, selected_group_keys,
            custom_vars_str="", include_mortality=False, output_path="",
            convert_units=False, auto_aggregate=False, profile=None,
            skip_cleanup=False):
        """
        执行完整流程
        - selected_cycle_suffixes: 周期后缀列表 ['H','I']
        - selected_group_keys: 变量组key列表 ['demo','bmx','smq']
        - custom_vars_str: 自定义变量字符串，格式 "BMI,DEMO:RIAGENDR" 或 "RIAGENDR,RIDAGEYR"
            如果没有指定表名，用所选组涉及的表；也可指定表名如 "DEMO:RIAGENDR"
        - include_mortality: 是否包含死亡率数据
        - output_path: 输出路径
        """
        self.is_running = True

        try:
            self.log("="*60)
            self.log(f"NHANES 数据处理开始")
            cycles_str = ", ".join([f"{c[0]}({c[2]})" for c in NHANES_CYCLES
                                     if c[2] in selected_cycle_suffixes])
            self.log(f"选择的周期: {cycles_str}")
            self.log(f"选择的变量组: {', '.join(selected_group_keys)}")
            if custom_vars_str.strip():
                self.log(f"自定义变量: {custom_vars_str.strip()}")
            if include_mortality:
                self.log("包含死亡率数据")
            self.log("="*60)

            # ---- 0. 依赖检查 ----
            if profile is not None:
                dep_warnings = check_dependencies(selected_cycle_suffixes, selected_group_keys, profile)
                for w in dep_warnings:
                    self.log(w)
                if dep_warnings:
                    self.log("⛔ 请修正上述配置问题后重试")
                    raise ValueError("\n".join(dep_warnings))

            # ---- 1. 收集所有需要下载的表（跨全部已选组去重）----
            all_needed_tables = set()
            for g in VARIABLE_GROUPS:
                if g["group_key"] in selected_group_keys:
                    all_needed_tables.add(g["table_name"])
            # 补充关联表
            lipid_related = {"TCHOL", "HDL", "TRIGLY"}
            for tn in list(all_needed_tables):
                if tn in lipid_related:
                    all_needed_tables.update(lipid_related)
            # 血常规组需额外下载HSCRP表（CRP所在表，E/F周期可用，G周期无）
            if "CBC" in all_needed_tables:
                all_needed_tables.add("HSCRP")
                self.log("  ℹ️  已添加HSCRP表（CRP数据，E/F周期可用）")
            all_needed_tables.discard("MORT")

            self.log(f"需要下载的表格: {sorted(all_needed_tables)}")

            # ---- 2. 下载XPT数据 ----
            all_raw = {}  # suffix -> {table_name -> df}
            # 先统计需要检查的文件
            total_needed = 0
            cached_count = 0
            to_download = []  # (suffix, cycle_info, tn, try_tn)
            for suffix in selected_cycle_suffixes:
                cycle_info = None
                for c in NHANES_CYCLES:
                    if c[2] == suffix:
                        cycle_info = c
                        break
                if not cycle_info:
                    continue
                cycle_name, start_year, _ = cycle_info
                all_raw[suffix] = {}
                for tn in sorted(all_needed_tables):
                    table_candidates = TABLE_ALIASES.get(tn, [tn])
                    for try_tn in table_candidates:
                        total_needed += 1
                        # 检查缓存
                        fname = f"{try_tn}_{suffix}.xpt" if suffix != "A" else f"{try_tn}.xpt"
                        cache_path = os.path.join(self.downloader.cache_dir, fname)
                        if os.path.exists(cache_path) and os.path.getsize(cache_path) > 1024:
                            cached_count += 1
                            # 缓存命中，直接读
                            df = read_xpt(cache_path)
                            if df is not None:
                                all_raw[suffix][tn] = df
                            break
                        else:
                            to_download.append((suffix, cycle_info, tn, try_tn))
            self.log(f"缓存: {cached_count}/{total_needed} 个文件已存在")
            if to_download:
                self.log(f"需要下载: {len(to_download)} 个文件")
                for suffix, cycle_info, tn, try_tn in to_download:
                    cycle_name, start_year, _ = cycle_info
                    # 修复2: 遍历所有别名，第一个成功即停止
                    candidates = TABLE_ALIASES.get(tn, [tn])
                    start_idx = candidates.index(try_tn) if try_tn in candidates else 0
                    downloaded_ok = False
                    for ci in range(start_idx, len(candidates)):
                        alt_tn = candidates[ci]
                        # 如果已缓存，直接读
                        fname = f"{alt_tn}_{suffix}.xpt" if suffix != "A" else f"{alt_tn}.xpt"
                        cache_path = os.path.join(self.downloader.cache_dir, fname)
                        if os.path.exists(cache_path) and os.path.getsize(cache_path) > 1024:
                            df = read_xpt(cache_path)
                            if df is not None:
                                all_raw[suffix][tn] = df
                                self.log(f"   ✅ {try_tn}_{suffix} 缓存回退到 {alt_tn}_{suffix} ({len(df)}行)")
                                downloaded_ok = True
                                break
                        # 下载
                        self.log(f"  下载: {alt_tn}_{suffix}")
                        fp, ok = self.downloader.download_xpt(alt_tn, suffix, start_year, self.log)
                        if ok and fp:
                            df = read_xpt(fp)
                            if df is not None:
                                all_raw[suffix][tn] = df
                                self.log(f"    {len(df)} 行, {len(df.columns)} 列")
                                downloaded_ok = True
                                break
                        else:
                            self.log(f"   ⚠️ {alt_tn}_{suffix} 下载失败，尝试下一别名")
                    if not downloaded_ok:
                        self.log(f"   ❌ {tn}_{suffix} 所有别名均下载失败（可能CDC无此表）")
            else:
                self.log("所有文件均已缓存，跳过下载")
                # 从缓存读取未加载的表
                for suffix, cycle_info, tn, try_tn in [(s, None, t, tt) 
                    for s in selected_cycle_suffixes
                    for t in sorted(all_needed_tables)
                    for tt in TABLE_ALIASES.get(t, [t])]:
                    if tn in all_raw.get(suffix, {}):
                        continue
                    fname = f"{tt}_{suffix}.xpt" if suffix != "A" else f"{tt}.xpt"
                    cache_path = os.path.join(self.downloader.cache_dir, fname)
                    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 1024:
                        df = read_xpt(cache_path)
                        if df is not None:
                            all_raw.setdefault(suffix, {})[tn] = df

            # ---- 3. 下载死亡率数据(如需) ----
            mort_data = {}  # suffix -> df
            if include_mortality:
                self.log(f"\n下载死亡率数据...")
                for suffix in selected_cycle_suffixes:
                    cycle_info = None
                    for c in NHANES_CYCLES:
                        if c[2] == suffix:
                            cycle_info = c
                            break
                    if not cycle_info:
                        continue
                    _, start_year, _ = cycle_info
                    fp, ok = self.downloader.download_mortality(suffix, start_year, self.log)
                    if ok and fp:
                        mdf = read_mortality(fp)
                        if mdf is not None:
                            mort_data[suffix] = mdf
                            self.log(f"    {cycle_info[0]}: {len(mdf)} 条死亡率记录")

            # ---- 3b. 收集原始RXQ_RX数据（用于药物排除）----
            rxq_raw_data = None
            if profile is not None and profile.config["exclusion"].get("medication_classes"):
                rxq_frames = []
                for suffix in selected_cycle_suffixes:
                    if suffix in all_raw and "RXQ_RX" in all_raw[suffix]:
                        rxq_frames.append(all_raw[suffix]["RXQ_RX"])
                if rxq_frames:
                    rxq_raw_data = pd.concat(rxq_frames, ignore_index=True)
                    self.log(f"收集 RXQ_RX 数据: {len(rxq_raw_data)} 条, "
                             f"涉及 {rxq_raw_data['SEQN'].nunique() if 'SEQN' in rxq_raw_data.columns else 0} 人")
                else:
                    self.log("⚠️ 需要药物排除但未找到 RXQ_RX 数据，请确保选择了'处方药'组")

            # ---- 4. 数据处理与合并 ----
            self.log(f"\n处理数据...")
            processed = {}
            per_cycle_var_log = {}  # suffix -> {col_label: {na_pct, n}}
            cycle_col_names = {}    # suffix -> [list of column labels]

            for suffix in selected_cycle_suffixes:
                cycle_info = None
                for c in NHANES_CYCLES:
                    if c[2] == suffix:
                        cycle_info = c
                        break
                if not cycle_info:
                    continue
                cycle_name = cycle_info[0]
                raw = all_raw.get(suffix, {})
                if not raw:
                    continue

                # 从DEMO表获取SEQN作为基础
                demo_df = raw.get("DEMO")
                if demo_df is None or "SEQN" not in demo_df.columns:
                    self.log(f"   ⚠️ {cycle_name} 缺少DEMO表，跳过")
                    continue

                base = demo_df[["SEQN"]].copy()
                base["序号(SEQN)"] = pd.to_numeric(base["SEQN"], errors='coerce').astype('Int64')
                base.drop(columns=["SEQN"], inplace=True)
                base["调查周期"] = cycle_name

                # 收集所有已下载表的列名索引 {col_name: {table_name, df}}
                all_table_index = {}
                for tn, tdf in raw.items():
                    for col in tdf.columns:
                        if col not in all_table_index:
                            all_table_index[col] = []
                        all_table_index[col].append((tn, tdf))

                # ---- 预聚合：RXQ_RX（1:N表，需先聚合成1条/人）----
                rxq_aggregated = None
                if "rxq" in selected_group_keys and "RXQ_RX" in raw:
                    rxq_raw = raw["RXQ_RX"]
                    if "SEQN" in rxq_raw.columns and "RXDUSE" in rxq_raw.columns:
                        rxq_agg = rxq_raw.groupby("SEQN", as_index=False).agg(
                            RXDUSE=("RXDUSE", "first"),
                            RXDCOUNT=("RXDCOUNT", "max"),
                        )
                        rxq_agg["使用处方药(Rx-Use)"] = rxq_agg["RXDUSE"].map({1: "是", 2: "否"})
                        rxq_agg["处方药数(Rx-Count,n)"] = rxq_agg["RXDCOUNT"]
                        rxq_aggregated = rxq_agg[["SEQN", "使用处方药(Rx-Use)", "处方药数(Rx-Count,n)"]].copy()
                        rxq_aggregated["序号(SEQN)"] = pd.to_numeric(rxq_aggregated["SEQN"], errors="coerce").astype("Int64")
                        self.log(f"   📊 RXQ_RX 已聚合: {len(rxq_raw)} 条→{len(rxq_aggregated)} 人")

                # 标记1:N表（不参与标准左连接，避免行膨胀）
                ONE_TO_MANY_TABLES = {"RXQ_RX"}

                # 按组处理（跨表搜索变量）
                for g in VARIABLE_GROUPS:
                    gk = g["group_key"]
                    if gk not in selected_group_keys:
                        continue

                    # 跳过1:N表的组（RXQ_RX已在上方预聚合）
                    if g.get("table_name", "") in ONE_TO_MANY_TABLES:
                        continue

                    group_vars = [v for v in g["variables"] if not v.get("is_id")]
                    if not group_vars:
                        continue

                    # 为本组每个变量找到数据
                    vars_found = {}  # canonical_name -> (series, var_info, mapping_or_None)
                    for v in group_vars:
                        vn = v["var_name"]
                        candidates = [vn] + VARIABLE_ALIASES.get(vn, [])
                        found_col = None
                        found_df = None
                        for cname in candidates:
                            if cname in all_table_index:
                                # 取第一个匹配的表
                                found_col = cname
                                found_df = all_table_index[cname][0][1]
                                break
                        if found_col and found_df is not None:
                            series = found_df[found_col].copy()
                            # 传入变量名和标签，以便对特殊变量（如睡眠时长）进行正确处理
                            var_label = v.get("var_label", "")
                            var_name = v.get("var_name", "")
                            series = apply_missing_clean(series, var_label=var_label, var_name=var_name)
                            if "mapping" in v:
                                series = series.map(v["mapping"])
                            vars_found[vn] = (series, v["var_label"], v.get("mapping"))

                    if not vars_found:
                        continue

                    # 合并到base
                    merge_df = pd.DataFrame({"序号(SEQN)": pd.to_numeric(demo_df["SEQN"], errors='coerce').astype('Int64')})
                    for vn, (series, label, _) in vars_found.items():
                        # 需要确保索引对齐
                        original_df = [t for t in all_table_index.values() if any(vn == x[0] or vn in VARIABLE_ALIASES.get(vn,[]) for x in t)]
                        # 简单方法：用series自身
                        s = series.copy()
                        merge_df[label] = s.values if len(s) == len(merge_df) else np.nan

                    # 通过SEQN精确左连接
                    # 需要从原始表获取带SEQN的数据框
                    merged_ok = False
                    for vn, (series, label, mapping) in vars_found.items():
                        # 找到这个变量来自哪个表
                        src_table = None
                        src_col = None
                        for cname in [vn] + VARIABLE_ALIASES.get(vn, []):
                            if cname in all_table_index:
                                src_table = all_table_index[cname][0][1]
                                src_col = cname
                                break
                        if src_table is not None and "SEQN" in src_table.columns:
                            # 修复1: 如果变量已合并过（如demo_core和demo_detail都取WTMEC2YR），跳过避免_x/_y后缀
                            if label in base.columns:
                                self.log(f"   ⏩ {cycle_name}: {gk} +{label} (已存在，跳过重复)")
                                continue
                            sub = src_table[["SEQN", src_col]].copy()
                            sub.rename(columns={"SEQN": "序号(SEQN)", src_col: label}, inplace=True)
                            sub["序号(SEQN)"] = pd.to_numeric(sub["序号(SEQN)"], errors='coerce').astype('Int64')
                            # 传入变量标签以便对特殊变量正确处理
                            sub[label] = apply_missing_clean(sub[label], var_label=label, var_name=src_col)
                            # 重新应用变量映射（如果有），确保值编码正确
                            if mapping:
                                sub[label] = sub[label].map(mapping)
                            before = len(base.columns)
                            base = base.merge(sub[["序号(SEQN)", label]], on="序号(SEQN)", how="left")
                            if len(base.columns) > before:
                                merged_ok = True
                                self.log(f"   ✅ {cycle_name}: {gk} +{label}")

                # 处理自定义变量（跨表搜索）
                custom_vars = [c.strip() for c in custom_vars_str.split(",") if c.strip()]
                for cv_name in custom_vars:
                    cv_upper = cv_name.upper()
                    found = False
                    for tn, tdf in raw.items():
                        for col in tdf.columns:
                            if col.upper() == cv_upper:
                                sub = tdf[["SEQN", col]].copy()
                                sub.rename(columns={"SEQN": "序号(SEQN)", col: cv_upper}, inplace=True)
                                sub["序号(SEQN)"] = pd.to_numeric(sub["序号(SEQN)"], errors='coerce').astype('Int64')
                                sub[cv_upper] = apply_missing_clean(sub[cv_upper])
                                before = len(base.columns)
                                base = base.merge(sub[["序号(SEQN)", cv_upper]], on="序号(SEQN)", how="left")
                                if len(base.columns) > before:
                                    self.log(f"   ✅ {cycle_name}: 自定义 +{cv_upper}")
                                found = True
                                break
                        if found:
                            break
                    if not found:
                        self.log(f"   ⏭️  {cycle_name}: 自定义变量 '{cv_name}' 未找到")

                # ---- 合并预聚合的RXQ_RX（1条/人）----
                if rxq_aggregated is not None:
                    merge_cols = [c for c in rxq_aggregated.columns if c not in ("SEQN", "RXDUSE", "RXDCOUNT", "序号(SEQN)")]
                    if len(merge_cols) >= 1:
                        before = len(base.columns)
                        base = base.merge(rxq_aggregated[["序号(SEQN)"] + merge_cols], on="序号(SEQN)", how="left")
                        added = len(base.columns) - before
                        self.log(f"   ✅ {cycle_name}: RXQ_RX(聚合) +{added} 变量")

                # 合并死亡率数据
                if include_mortality and suffix in mort_data:
                    mdf = mort_data[suffix]
                    mdf_subset = mdf[["SEQN"]].copy()
                    mdf_subset["序号(SEQN)"] = pd.to_numeric(mdf_subset["SEQN"], errors='coerce').astype('Int64')
                    mdf_subset.drop(columns=["SEQN"], inplace=True)

                    for v in MORTALITY_GROUPS[0]["variables"]:
                        vn = v["var_name"]
                        if vn in mdf.columns:
                            col = mdf[vn].copy()
                            col = apply_missing_clean(col)
                            if "mapping" in v:
                                col = col.map(v["mapping"])
                            mdf_subset[v["var_label"]] = col

                    merge_col = "序号(SEQN)"
                    right_cols = [c for c in mdf_subset.columns if c != merge_col]
                    if right_cols:
                        before = len(base.columns)
                        base = base.merge(mdf_subset[[merge_col] + right_cols], on=merge_col, how="left")
                        added = len(base.columns) - before
                        self.log(f"   ✅ {cycle_name}: 死亡率 +{added} 变量")

                processed[suffix] = base
                self.log(f"   📊 {cycle_name}: {len(base)} 条, {len(base.columns)} 列")
                # 记录该周期变量存在性
                cycle_log = {}
                for c in base.columns:
                    if c in ("序号(SEQN)", "调查周期"):
                        continue
                    na_pct = base[c].isna().mean() * 100
                    cycle_log[c] = {"na_pct": round(na_pct, 1), "n": int(base[c].notna().sum())}
                per_cycle_var_log[suffix] = cycle_log
                cycle_col_names[suffix] = list(base.columns)

            # ---- 5. 合并所有周期 ----
            self.log(f"\n合并周期...")
            if not processed:
                raise Exception("没有成功处理任何数据！")
            final = pd.concat(processed.values(), ignore_index=True)
            self.log(f"合并完成: {len(final)} 条, {len(final.columns)} 列")

            # ---- 5b. 单位换算 ----
            if convert_units:
                self.log(f"\n应用单位换算(美制→中国常用单位)...")
                before_cols = len(final.columns)
                final = apply_unit_conversion(final, enabled=True)
                added = len(final.columns) - before_cols
                if added > 0:
                    self.log(f"✅ 新增 {added} 列换算后指标")

            # ---- 自动聚合细碎变量 ----
            if auto_aggregate:
                self.log(f"\n自动聚合细碎变量(合并为论文可用指标)...")
                agg = DataAggregator(final)
                final = agg.run()
                for msg in agg.logs:
                    self.log(f"  {msg}")

            # ---- 清洗管道（纳入排除+衍生变量+异常值标记）----
            reporter = LogReporter()
            if profile is not None:
                self.log(f"\n应用研究方案清洗: {profile.config.get('study_name','')}")
                pipeline = CleaningPipeline(final, profile, rxq_data=rxq_raw_data, reporter=reporter)
                final = pipeline.run()
                for msg in pipeline.logs:
                    self.log(f"  {msg}")
                if pipeline.removed_count > 0:
                    self.log(f"  清洗后剩余: {len(final)} 人")

            # ---- 质控报告 ----
            qc_report = {}
            qc_report["总记录数"] = len(final)
            qc_report["总变量数"] = len(final.columns)
            qc_report["周期分布"] = {}
            if "调查周期" in final.columns:
                qc_report["周期分布"] = final["调查周期"].value_counts().to_dict()
            qc_report["变量缺失率"] = {}
            for col in final.columns:
                if col in ("序号(SEQN)", "调查周期"):
                    continue
                missing_pct = final[col].isna().mean() * 100
                if missing_pct > 0:
                    qc_report["变量缺失率"][col] = round(missing_pct, 1)
            # 数值变量范围检查
            qc_report["数值范围"] = {}
            qc_report["偏态变量"] = {}  # 新增：偏态变量标记
            numeric_cols = final.select_dtypes(include=[np.number]).columns
            for col in numeric_cols:
                if col in ("序号(SEQN)",):
                    continue
                valid = final[col].dropna()
                if len(valid) > 0:
                    # 基础统计
                    qc_report["数值范围"][col] = {
                        "最小值": round(float(valid.min()), 2),
                        "最大值": round(float(valid.max()), 2),
                        "均值": round(float(valid.mean()), 2),
                        "标准差": round(float(valid.std()), 2),
                        "中位数": round(float(valid.median()), 2),
                        "Q1(25%)": round(float(np.percentile(valid, 25)), 2),
                        "Q3(75%)": round(float(np.percentile(valid, 75)), 2),
                        "IQR": round(float(np.percentile(valid, 75) - np.percentile(valid, 25)), 2),
                    }
                    
                    # 偏态检测（偏度>1或<-1视为偏态）
                    try:
                        skewness = valid.skew()
                        qc_report["数值范围"][col]["偏度"] = round(float(skewness), 2)
                        
                        if abs(skewness) > 1:
                            # 标记为偏态变量
                            qc_report["偏态变量"][col] = {
                                "偏度": round(float(skewness), 2),
                                "类型": "右偏态" if skewness > 1 else "左偏态",
                                "建议": "报告中位数(IQR)而非均值±标准差；参数检验前需对数转换",
                                "中位数": round(float(valid.median()), 2),
                                "IQR": f"{round(float(np.percentile(valid, 25)), 2)}-{round(float(np.percentile(valid, 75)), 2)}",
                            }
                    except:
                        pass

            # 排序
            if "序号(SEQN)" in final.columns:
                final = final.sort_values("序号(SEQN)").reset_index(drop=True)

            # ---- 列清理（DB模式跳过，保留全部列）----
            if not skip_cleanup:
                # 构建动态白名单：用户所选变量组中所有变量的标签
                keep_labels = set()
                for g in VARIABLE_GROUPS:
                    if g["group_key"] in selected_group_keys:
                        for v in g["variables"]:
                            if not v.get("is_id"):
                                keep_labels.add(v["var_label"])
                # 补充死亡率变量标签（如果启用了死亡率）
                if include_mortality and "mort" in selected_group_keys:
                    for g in MORTALITY_GROUPS:
                        for v in g["variables"]:
                            keep_labels.add(v["var_label"])
                final = clean_and_reorder_columns(final, keep_labels=keep_labels)
                self.log(f"\n列清理后: {len(final.columns)} 列（动态白名单: {len(keep_labels)}个变量标签）")

            # ---- 变量元数据生成 ----
            if output_path:
                import json
                metadata_rows = []
                for c in final.columns:
                    if c in ("序号(SEQN)", "调查周期"):
                        continue
                    row = {"变量名": c}
                    for suffix in selected_cycle_suffixes:
                        cycle_name = [ci[0] for ci in NHANES_CYCLES if ci[2] == suffix][0]
                        if suffix in per_cycle_var_log and c in per_cycle_var_log[suffix]:
                            info = per_cycle_var_log[suffix][c]
                            row[f"{cycle_name}_缺失率%"] = info["na_pct"]
                            row[f"{cycle_name}_有效N"] = info["n"]
                        else:
                            row[f"{cycle_name}_缺失率%"] = -1
                            row[f"{cycle_name}_有效N"] = 0
                    metadata_rows.append(row)
                md_df = pd.DataFrame(metadata_rows)
                # 添加短名
                short_names = {}
                for c in final.columns:
                    if c in self.COLUMN_MAP_963:
                        short_names[c] = self.COLUMN_MAP_963[c]
                    else:
                        m = re.search(r'\(([A-Za-z][A-Za-z0-9_\-]+)', c)
                        short_names[c] = m.group(1) if m else c
                md_df["短名"] = md_df["变量名"].map(short_names)
                # 排序：高缺失率优先
                na_cols = [f"{ci[0]}_缺失率%" for ci in NHANES_CYCLES if ci[2] in selected_cycle_suffixes]
                if na_cols:
                    md_df = md_df.sort_values(by=na_cols, ascending=False)
                self._cycle_var_metadata = md_df  # 暂存供后续使用
            
            # ---- 6. 导出 ----
            full_csv = None
            if output_path:
                # 计算年份范围作为文件夹名
                cycle_names = []
                for c in NHANES_CYCLES:
                    if c[2] in selected_cycle_suffixes:
                        cycle_names.append(c[0])
                if cycle_names:
                    folder_name = f"{cycle_names[0].split('-')[0]}-{cycle_names[-1].split('-')[1]}"
                else:
                    folder_name = "unknown"

                if os.path.isdir(output_path):
                    # 用户只选了文件夹 → 自动生成文件名
                    base_folder = output_path
                    base_name = "NHANES_" + "_".join(selected_cycle_suffixes)
                else:
                    # 兼容旧版：从文件路径提取
                    base_folder = os.path.dirname(os.path.abspath(output_path))
                    base_name = os.path.splitext(os.path.basename(output_path))[0]
                    base_name = re.sub(r'_\d{4}-\d{4}$', '', base_name)
                
                out_dir = os.path.join(base_folder, folder_name)
                os.makedirs(out_dir, exist_ok=True)
                
                # ---- 仅1份文件：未加权数据（短列名 + 全列保留）----
                df_short = self._make_short_names(final)
                unweighted_csv = os.path.join(out_dir, f"{base_name}_未加权数据.csv")
                df_short.to_csv(unweighted_csv, index=False, encoding="utf-8-sig")
                self.log(f"  ✅ 未加权数据: {unweighted_csv} ({len(df_short)}人, {len(df_short.columns)}列)")
                self.log(f"     ⚠️ 加权分析请使用R survey包: svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR)")

                # ---- 变量元数据与数据完整性报告 ----
                if hasattr(self, '_cycle_var_metadata') and self._cycle_var_metadata is not None:
                    md_path = os.path.join(out_dir, f"{base_name}_变量元数据.csv")
                    self._cycle_var_metadata.to_csv(md_path, index=False, encoding="utf-8-sig")
                    self.log(f"  📋 变量元数据: {md_path}")
                    
                    # 数据完整性报告
                    rep_lines = ["="*60]
                    rep_lines.append("NHANES 数据完整性报告")
                    rep_lines.append(f"周期: {', '.join(selected_cycle_suffixes)}")
                    rep_lines.append("="*60)
                    rep_lines.append("")
                    for suffix in selected_cycle_suffixes:
                        cycle_name = [ci[0] for ci in NHANES_CYCLES if ci[2] == suffix][0]
                        rep_lines.append(f"--- {cycle_name} ---")
                        if suffix not in per_cycle_var_log:
                            rep_lines.append(f"  ❌ 无数据"); continue
                        log = per_cycle_var_log[suffix]
                        # 按缺失率分组
                        low_na = [(k,v) for k,v in log.items() if v['na_pct'] <= 30]
                        high_na = [(k,v) for k,v in log.items() if v['na_pct'] > 70]
                        mid_na = [(k,v) for k,v in log.items() if 30 < v['na_pct'] <= 70]
                        rep_lines.append(f"  总变量: {len(log)}")
                        rep_lines.append(f"  低缺失(≤30%): {len(low_na)}")
                        if high_na:
                            rep_lines.append(f"  高缺失(>70%): {len(high_na)}")
                            for k,v in sorted(high_na, key=lambda x:-x[1]['na_pct'])[:5]:
                                rep_lines.append(f"    - {k}: {v['na_pct']}% (有效N={v['n']})")
                        if mid_na:
                            rep_lines.append(f"  中度缺失(30-70%): {len(mid_na)}")
                        rep_lines.append("")
                    # 跨周期差异
                    if len(selected_cycle_suffixes) >= 2:
                        rep_lines.append("--- 跨周期变量差异 ---")
                        ref_suf = selected_cycle_suffixes[0]
                        ref_vars = set(cycle_col_names.get(ref_suf, []))
                        for suf in selected_cycle_suffixes[1:]:
                            cur_vars = set(cycle_col_names.get(suf, []))
                            missing = ref_vars - cur_vars
                            extra = cur_vars - ref_vars
                            cyc_name = [ci[0] for ci in NHANES_CYCLES if ci[2] == suf][0]
                            if missing:
                                rep_lines.append(f"  {cyc_name} 缺失（vs {ref_suf}）: {', '.join(sorted(missing)[:10])}")
                            if extra:
                                rep_lines.append(f"  {cyc_name} 新增: {', '.join(sorted(extra)[:10])}")
                        rep_lines.append("")
                    rep_lines.append("--- 缺失原因分类 ---")
                    rep_lines.append("1. NA = 此变量在该周期被测量但个体无数据（HDL仅检测空腹等）")
                    rep_lines.append("2. NA%=-1 = 此变量在该周期无对应检测（如CRP在G周期不存在）")
                    rep_lines.append("3. 跨周期差异 = CDC在不同年份增加了/删除了检测项目")
                    rep_lines.append("")
                    integrity_path = os.path.join(out_dir, f"{base_name}_数据完整性报告.txt")
                    with open(integrity_path, "w", encoding="utf-8") as f:
                        f.write("\n".join(rep_lines))
                    self.log(f"  📄 数据完整性报告: {integrity_path}")
                
                # 筛选报告（详细步骤）
                if profile is not None:
                    cycle_dist_str = ", ".join([f"{k}={v}" for k,v in qc_report.get("周期分布",{}).items()]) if qc_report.get("周期分布") else ""
                    reporter.qc_report = qc_report
                    report_path = reporter.write_txt(
                        unweighted_csv,
                        len(final), len(final.columns), cycle_dist_str
                    )
                    if report_path:
                        self.log(f"  📄 筛选报告: {report_path}")
                
                full_csv = unweighted_csv  # 返回值兼容
            else:
                self.log(f"\n(跳过导出，仅返回DataFrame)")

            self.log(f"\n✅ 完成! 共 {len(final)} 条记录, {len(final.columns)} 列")

            # 打印质控摘要
            self.log(f"\n📊 质控摘要:")
            self.log(f"   周期: {', '.join([f'{k}={v}' for k,v in qc_report.get('周期分布',{}).items()])}")
            high_missing = {k:v for k,v in qc_report.get("变量缺失率",{}).items() if v > 50}
            if high_missing:
                self.log(f"   ⚠️ 高缺失率变量(>50%): {len(high_missing)} 个")
                for k,v in list(high_missing.items())[:5]:
                    self.log(f"     - {k}: {v}%")
            self.log(f"   建议打开「预览结果」查看完整质控报告")

            # ---- 自动数据健康检查 ----
            try:
                from qc_engine import QCEngine
                qc = QCEngine(final)
                qc_report_raw = qc.run_all()
                qc_passed = qc_report_raw['passed']
                qc_summary = qc_report_raw['summary']
                self.log(f"   {'✅' if qc_passed else '❌'} 数据健康检查: {qc_summary}")
                if not qc_passed:
                    fails = [c for c in qc_report_raw['checks'] if not c['passed']]
                    for f in fails[:3]:
                        self.log(f"      ❌ {f['name']}: {f['detail'][:60]}")
            except Exception as e:
                self.log(f"   ⚠️ 健康检查跳过: {e}")

            return {"success": True, "rows": len(final), "cols": len(final.columns),
                    "file_path": full_csv if output_path else None, "df": final, "qc_report": qc_report}

        except Exception as e:
            import traceback
            self.log(f"\n❌ 错误: {e}")
            self.log(traceback.format_exc())
            return {"success": False, "error": str(e)}
        finally:
            self.is_running = False

    def run_async(self, *args, **kwargs):
        t = threading.Thread(target=self.run, args=args, kwargs=kwargs, daemon=True)
        t.start()
        return t

    def _export_weighted_to(self, df, target_path, seed=42):
        """
        保存带有权重列的完整数据集（不含行复制）
        CDC规范：权重列(WTMEC2YR等)代表抽样系数，不应复制行来"模拟加权"
        正确的加权分析应在R的survey包或Python的statsmodels中完成
        """
        try:
            # 查找权重列（扩展匹配）
            wt_col = None
            all_cols = list(df.columns)
            candidates = ["WTSA2YR", "甲状腺权重(WT-Subset)", "甲状腺权重(WT-Subset)_x", "甲状腺权重(WT-Subset)_y", 
                          "体检权重(WT-MEC)", "WTMEC2YR",
                          "访谈权重(WT-Int)", "WTINT2YR", "wt_mec", "wt_int"]
            for c in candidates:
                if c in df.columns:
                    wt_col = c
                    break
            # 模糊匹配
            if not wt_col:
                for c in all_cols:
                    cu = c.upper().replace(" ","").replace("-","").replace("(","").replace(")","")
                    if "权重" in c or "WTMEC" in cu or "WTINT" in cu:
                        wt_col = c
                        break
            if not wt_col:
                self.log(f"   ⏭️ 加权版跳过: 未找到权重列")
                return None

            # ✅ 正确做法：直接保存原始数据，保留权重列
            # 不再使用 np.repeat 行复制（那是致命错误——虚假膨胀样本量）
            df.to_csv(target_path, index=False, encoding="utf-8-sig")
            self.log(f"  ✅ 加权版(保留权重列): {target_path} ({len(df)}人, {len(df.columns)}列, 权重列: {wt_col})")
            self.log(f"     ⚠️ 加权分析请使用R survey包: svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~{wt_col})")
            return target_path
        except Exception as e:
            self.log(f"  ⚠️ 加权版保存失败: {e}")
            return None

    # ========== 9.63分析工具兼容版 ==========
    COLUMN_MAP_963 = {
        "序号(SEQN)": "序号(SEQN)",
        "调查周期": "调查周期",
        "性别(Gender)": "性别",
        "年龄(Age,岁)": "年龄",
        "BMI(kg/m²)": "BMI",
        "腰围(WC,cm)": "WC",
        "收缩压-1(SBP1,mmHg)": "SBP1",
        "舒张压-1(DBP1,mmHg)": "DBP1",
        "PSU(SDMVPSU)": "SDMVPSU",
        "分层(SDMVSTRA)": "SDMVSTRA",
        "体检权重(WT-MEC)": "WT-MEC",
        "访谈权重(WT-Int)": "WT-Int",
        "促甲状腺激素(TSH,mIU/L)": "TSH",
        "游离三碘甲腺原氨酸(FT3,pmol/L)": "FT3",
        "游离甲状腺素(FT4,pmol/L)": "FT4",
        "总甲状腺素(TT4,nmol/L)": "TT4",
        "总三碘甲腺原氨酸(TT3,nmol/L)": "TT3",
        "总胆固醇(TC,mmol/L)": "TC",
        "甘油三酯(TG,mmol/L)": "TG",
        "高密度脂蛋白(HDL,mmol/L)": "HDL-C",
        "低密度脂蛋白(LDL,mmol/L)": "LDL-C",
        "空腹血糖(Glu,mmol/L)": "Glu",
        "谷丙转氨酶(ALT,IU/L)": "ALT",
        "谷草转氨酶(AST,IU/L)": "AST",
        "γ-谷氨酰转移酶(GGT,IU/L)": "GGT",
        "肌酐(Cr,μmol/L)": "Cr",
        "尿酸(UA,μmol/L)": "UA",
        "尿素氮(BUN,mmol/L)": "BUN",
        "脉率(PR,bpm)": "PR",
        "白细胞(WBC,×10⁹/L)": "WBC",
        "血小板(PLT,×10⁹/L)": "PLT",
        "血红蛋白(Hb,g/L)": "Hb",
        "白蛋白(ALB,g/L)": "ALB",
        "总胆红素(TBIL,μmol/L)": "TBIL",
    }

    def _make_963_df(self, df):
        """
        将完整列名DataFrame转换为9.63兼容的短列名版本
        """
        rename_map = {}
        keep_cols = []
        for col in df.columns:
            if col in self.COLUMN_MAP_963:
                rename_map[col] = self.COLUMN_MAP_963[col]
                keep_cols.append(col)
            elif col in ("调查周期", "年份", "序号(SEQN)"):
                keep_cols.append(col)

        df963 = df[keep_cols].copy()
        df963.rename(columns=rename_map, inplace=True)

        # 确保年份列存在
        if "年份" not in df963.columns and "调查周期" in df.columns:
            year_map = {"2007-2008": 2007, "2009-2010": 2009, "2011-2012": 2011,
                        "2013-2014": 2013, "2015-2016": 2015, "2017-2018": 2017,
                        "2019-2020": 2019, "2021-2022": 2021, "2023-2024": 2023}
            df963["年份"] = df["调查周期"].map(year_map)

        # 性别映射：NHANES用1=男/2=女 → 9.63要求1=男/0=女
        if "性别" in df963.columns:
            s = pd.to_numeric(df963["性别"], errors="coerce")
            # 如果存在2（女），才做映射
            if s.eq(2).any():
                df963["性别"] = s.map({1: 1, 2: 0}).fillna(s)
            else:
                # 已是0/1格式，确保整数
                df963["性别"] = (s == 1).astype(int)

        return df963

    def _make_short_names(self, df):
        """
        将完整列名统一转为短名格式（保留所有列）
        '促甲状腺激素(TSH,mIU/L)' → 'TSH'
        '收缩压-1(SBP1,mmHg)' → 'SBP1'
        '性别(Gender)' → '性别'（保持原值"男"/"女"）
        无法提取短名的列保留原名
        """
        rename_map = {}
        for col in df.columns:
            # 已有映射的优先使用
            if col in self.COLUMN_MAP_963:
                rename_map[col] = self.COLUMN_MAP_963[col]
                continue
            # 通用模式：提取括号中的首个英文缩写
            # 匹配 中文(CODE,unit) 或 中文(CODE) 中的 CODE
            m = re.search(r'\(([A-Za-z][A-Za-z0-9_\-]+)', col)
            if m:
                rename_map[col] = m.group(1)
            # 保底：列名包含"序号(SEQN)"类 → 提取SEQN
            m2 = re.search(r'\(([A-Za-z0-9_\-]+)\)', col)
            if m2 and col not in rename_map:
                rename_map[col] = m2.group(1)

        df_out = df.rename(columns=rename_map)
        return df_out


# ============================================================================
# ===== GUI 界面 =====
# ============================================================================

def _bind_mousewheel(canvas):
    """绑定鼠标滚轮到 Canvas"""
    def _on_mousewheel(event):
        canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
    def _on_mousewheel_linux(event):
        canvas.yview_scroll(-1 if event.num == 4 else 1, "units")
    canvas.bind("<MouseWheel>", _on_mousewheel)
    canvas.bind("<Button-4>", _on_mousewheel_linux, add="+")
    canvas.bind("<Button-5>", _on_mousewheel_linux, add="+")

# ============================================================================
# 数据自动聚合模块（将细碎变量自动合并为论文可用指标）
# ============================================================================
# NHANES 很多变量是散碎的单项（如每种鱼吃没吃），
# 此模块自动识别变量模式，聚合为有意义的分类变量。
# 用户在 GUI 勾选"自动聚合"即可无需手动处理。

AGGREGATION_RULES = {
    "鱼类摄入(Fish-Eat)": {
        "description": "过去30天是否吃鱼(汇总DRD350A-K)",
        "type": "any_yes",
        "source_pattern": r"^DRD350[A-K]$",
        "value_mapping": {1: "是", 2: "否"},
    },
    "鱼类摄入次数(Fish-Freq)": {
        "description": "过去30天吃鱼总次数",
        "type": "sum",
        "source_pattern": r"^DRD350[A-K]Q$",
    },
    "贝类摄入(Shellfish-Eat)": {
        "description": "过去30天是否吃贝类(汇总DRD350H-K)",
        "type": "any_yes",
        "source_pattern": r"^DRD350[H-K]$",
        "value_mapping": {1: "是", 2: "否"},
    },
    "具体鱼类品种数(Fish-Species)": {
        "description": "过去30天吃了几种鱼(DRD370A-R)",
        "type": "count_yes",
        "source_pattern": r"^DRD370[A-R]$",
        "value_mapping": {1: 1, 2: 0},
    },
    "高汞鱼类摄入(HighHg-Fish)": {
        "description": "是否吃过高汞鱼类(鲭鱼/鲨鱼/剑鱼等)",
        "type": "any_yes",
        "source_pattern": r"^DRD370[DMS]$",  # 示例
        "value_mapping": {1: "是", 2: "否"},
    },
    "吸烟状态(Smoke-Status)": {
        "description": "综合吸烟状态(从不/曾吸/目前吸)",
        "type": "derive_smoke",
        "variables_needed": ["SMQ020", "SMQ040"],
        "value_mapping": {},
    },
    "饮酒状态(Drink-Status)": {
        "description": "综合饮酒状态(从不/曾饮/目前饮)",
        "type": "derive_drink",
        "variables_needed": ["ALQ101", "ALQ110"],
        "value_mapping": {},
    },
}

class DataAggregator:
    """自动识别并聚合NHANES细碎变量"""

    def __init__(self, df, cycle=""):
        self.df = df
        self.cycle = cycle
        self.columns = set(df.columns)
        self.logs = []

    def log(self, msg):
        self.logs.append(msg)
        print(f"  [聚合] {msg}")

    def run(self):
        """对DataFrame执行所有可用的聚合"""
        before = len(self.df.columns)
        for rule_name, rule in AGGREGATION_RULES.items():
            try:
                self._apply_rule(rule_name, rule)
            except Exception as e:
                self.log(f"⚠️ {rule_name} 聚合失败: {e}")
        added = len(self.df.columns) - before
        if added > 0:
            self.log(f"✅ 共新增 {added} 个聚合变量")
        return self.df

    def _apply_rule(self, name, rule):
        rtype = rule["type"]

        if rtype == "any_yes":
            # 任一变量为"是"则结果为"是"
            pattern = rule["source_pattern"]
            matched = [c for c in self.columns if self._match(c, pattern)]
            if not matched:
                return
            # 检查值是否在列中
            result = pd.Series(False, index=self.df.index)
            for col in matched:
                if col in self.df.columns:
                    for yes_val in [1, "1"]:
                        result = result | (self.df[col] == yes_val)
            self.df[name] = result.map({True: "是", False: "否"})
            self.log(f"✅ {name}: 合并 {len(matched)} 个变量")

        elif rtype == "sum":
            pattern = rule["source_pattern"]
            matched = [c for c in self.columns if self._match(c, pattern)]
            if not matched:
                return
            total = pd.Series(0, dtype=float, index=self.df.index)
            valid_count = 0
            for col in matched:
                if col in self.df.columns:
                    numeric = pd.to_numeric(self.df[col], errors="coerce").fillna(0)
                    total = total + numeric
                    valid_count += 1
            if valid_count > 0:
                self.df[name] = total
                self.log(f"✅ {name}: 累计 {valid_count} 个数量变量")

        elif rtype == "count_yes":
            pattern = rule["source_pattern"]
            matched = [c for c in self.columns if self._match(c, pattern)]
            if not matched:
                return
            count = pd.Series(0, dtype=int, index=self.df.index)
            for col in matched:
                if col in self.df.columns:
                    count = count + (self.df[col] == 1).astype(int)
            self.df[name] = count
            self.log(f"✅ {name}: 汇总 {len(matched)} 个项目")

        elif rtype == "derive_smoke":
            needed = rule.get("variables_needed", [])
            # 判断是否有所需变量（可能已被翻译成中文标签）
            smq020 = self._find_col(needed[0]) or self._find_col("曾吸100支烟")
            smq040 = self._find_col(needed[1]) or self._find_col("目前吸烟频率")
            if smq020 is None or smq040 is None:
                return
            smoke = pd.Series("未知", index=self.df.index)
            never = (self.df[smq020] == 2) | (self.df[smq020] == "否")
            ever = (self.df[smq020] == 1) | (self.df[smq020] == "是")
            now = (self.df[smq040] == 1) | (self.df[smq040] == "每天吸")
            some = (self.df[smq040] == 2) | (self.df[smq040] == "偶尔吸")
            quit_ = (self.df[smq040] == 3) | (self.df[smq040] == "不吸")
            smoke.loc[never] = "从不吸烟"
            smoke.loc[ever & now] = "目前每天吸"
            smoke.loc[ever & some] = "目前偶尔吸"
            smoke.loc[ever & quit_] = "曾吸已戒"
            self.df[name] = smoke
            self.log(f"✅ {name}: 综合吸烟分类")

        elif rtype == "derive_drink":
            needed = rule.get("variables_needed", [])
            alq101 = self._find_col(needed[0]) or self._find_col("曾饮酒12次以上")
            alq110 = self._find_col(needed[1]) or self._find_col("过去1年饮酒")
            if alq101 is None:
                return
            drink = pd.Series("未知", index=self.df.index)
            never = (self.df[alq101] == 2) | (self.df[alq101] == "否")
            ever = (self.df[alq101] == 1) | (self.df[alq101] == "是")
            current = False
            if alq110:
                current = (self.df[alq110] == 1) | (self.df[alq110] == "是")
            drink.loc[never] = "从不饮酒"
            if alq110:
                drink.loc[ever & current] = "目前饮酒"
                drink.loc[ever & ~current] = "曾饮已戒"
            else:
                drink.loc[ever] = "曾饮酒"
            self.df[name] = drink
            self.log(f"✅ {name}: 综合饮酒分类")

    @staticmethod
    def _match(colname, pattern):
        """检查列名是否匹配正则模式"""
        import re
        return bool(re.match(pattern, colname.upper()))

    def _find_col(self, keyword):
        """按关键词在列名中查找，返回完整列名"""
        kw = keyword.upper()
        for c in self.columns:
            if kw in c.upper():
                return c
        return None


# ============================================================================
# 本地数据库构建器（一次性下载全部周期，保存为SQLite）
# ============================================================================

DB_DEFAULT_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "nhanes_local.db"
)

class LocalDBBuilder:
    """一键构建本地NHANES数据库：下载全部13个周期、全部变量组，存储为SQLite"""

    def __init__(self, db_path=DB_DEFAULT_PATH):
        self.db_path = db_path
        self.engine = NhanesEngine()
        self.log_callback = None
        self.is_running = False

    def log(self, msg):
        if self.log_callback:
            self.log_callback(msg)
        print(msg)

    def build(self, include_mortality=False, convert_units=True):
        """
        构建本地数据库
        1. 下载所有周期 + 所有表
        2. 提取所有变量组
        3. 标准化变量名
        4. 存入SQLite
        """
        self.is_running = True
        try:
            self.log("="*60)
            self.log("🚀 开始构建本地NHANES数据库")
            self.log(f"   目标: {self.db_path}")
            self.log(f"   周期: 全部13个 (1999-2024)")
            all_group_keys = [g["group_key"] for g in VARIABLE_GROUPS]
            self.log(f"   变量组: {len(all_group_keys)} 组, "
                     f"{sum(len(g['variables']) for g in VARIABLE_GROUPS)} 个变量")
            self.log(f"   死亡率: {'是' if include_mortality else '否'}")
            self.log(f"   单位换算: {'是' if convert_units else '否'}")
            self.log("="*60)

            all_suffixes = [c[2] for c in NHANES_CYCLES]
            result = self.engine.run(
                all_suffixes, all_group_keys,
                custom_vars_str="",
                include_mortality=include_mortality,
                output_path="",
                convert_units=False,       # DB存原始单位（美制），导出时再换算
                auto_aggregate=False,      # DB存原始变量，不聚合
                skip_cleanup=True)         # DB保留全部列，不删减

            if not result.get("success"):
                raise Exception(result.get("error", "处理失败"))

            df = result["df"]
            # 不应用 _make_short_names：DB保留原始长列名
            # 短列名转换在 query() 导出时再做
            
            self.log(f"\n💾 正在写入SQLite数据库...")

            import sqlite3
            conn = sqlite3.connect(self.db_path)
            df.to_sql("nhanes_all", conn, if_exists="replace", index=False)

            # 创建索引加速查询
            conn.execute("CREATE INDEX IF NOT EXISTS idx_seqn ON nhanes_all(\"序号(SEQN)\")")
            conn.execute("CREATE INDEX IF NOT EXISTS idx_cycle ON nhanes_all(\"调查周期\")")
            conn.commit()
            conn.close()

            self.log(f"✅ 数据库构建完成!")
            self.log(f"   路径: {self.db_path}")
            self.log(f"   记录: {len(df)} 条")
            self.log(f"   变量: {len(df.columns)} 列")
            self.log(f"   文件大小: {os.path.getsize(self.db_path)//1024//1024} MB")

            return {"success": True, "path": self.db_path, "rows": len(df), "cols": len(df.columns)}

        except Exception as e:
            import traceback
            self.log(f"\n❌ 构建失败: {e}")
            self.log(traceback.format_exc())
            return {"success": False, "error": str(e)}
        finally:
            self.is_running = False

    def query(self, cycles=None, variables=None, output_path="",
              convert_units=True, use_short_names=True):
        """从本地数据库查询数据
        
        Args:
            cycles: 周期列表，如 ['2007-2008', '2009-2010']，None=全部
            variables: 变量名列表，如 ['性别(Gender)', '年龄(Age,岁)']
            output_path: 导出文件路径
            convert_units: 是否做单位换算（美制→SI）
            use_short_names: 是否使用短列名（TSH/Glu/ALT等）
        """
        import sqlite3
        if not os.path.exists(self.db_path):
            return {"success": False, "error": f"数据库文件不存在: {self.db_path}"}

        conn = sqlite3.connect(self.db_path)
        query = 'SELECT * FROM nhanes_all'
        conditions = []
        if cycles:
            # DB里的周期列是中文格式 "2007-2008"，支持部分匹配
            cycle_conds = []
            for c in cycles:
                c = c.strip()
                if "-" in c:
                    cycle_conds.append(f"\"调查周期\" = '{c}'")
                else:
                    # 支持模糊匹配：输入 "2007" 匹配 "2007-2008"
                    cycle_conds.append(f"\"调查周期\" LIKE '{c}%'")
            if cycle_conds:
                conditions.append("(" + " OR ".join(cycle_conds) + ")")
        if conditions:
            query += " WHERE " + " AND ".join(conditions)

        try:
            df = pd.read_sql_query(query, conn)
            conn.close()

            # ---- 应用列清理（去除美标列+非核心列）----
            df = clean_and_reorder_columns(df)

            # ---- 单位换算（美制→中国常用单位）----
            if convert_units:
                df = apply_unit_conversion(df, enabled=True)

            # ---- 短列名转换 ----
            if use_short_names:
                df = self.engine._make_short_names(df)

            # ---- 用户选择的变量子集 ----
            if variables:
                available = [v for v in variables if v in df.columns]
                if "序号(SEQN)" in df.columns and "序号(SEQN)" not in available:
                    available.insert(0, "序号(SEQN)")
                if "调查周期" in df.columns and "调查周期" not in available:
                    available.insert(0, "调查周期")
                if available:
                    df = df[available]

            if output_path:
                if output_path.endswith(".xlsx"):
                    df.to_excel(output_path, index=False, engine="openpyxl")
                else:
                    df.to_csv(output_path, index=False, encoding="utf-8-sig")

            return {"success": True, "df": df, "rows": len(df), "cols": len(df.columns)}
        except Exception as e:
            if 'conn' in dir(): conn.close()
            return {"success": False, "error": str(e)}


# ============================================================================
# 药物分类系统（用于排除用药、筛选用药人群）
# ============================================================================
# RXQ_RX 表中 RXDDRUG 列包含药物名称，用关键词匹配分类
DRUG_CLASSIFICATION = {
    "降脂药(Lipid-lowering)": {
        "keywords": [
            "atorvastatin", "simvastatin", "rosuvastatin", "pravastatin",
            "lovastatin", "fluvastatin", "pitavastatin", "ezetimibe",
            "fenofibrate", "gemfibrozil", "niacin", "cholestyramine",
            "colesevelam", "colestipol", "evolocumab", "alirocumab",
            "icosapent", "omega-3", "vascepa", "lomitapide", "bempedoic",
            "阿托伐他汀", "辛伐他汀", "瑞舒伐他汀", "普伐他汀", "氟伐他汀", "匹伐他汀",
            "非诺贝特", "依折麦布", "益适纯", "可定", "立普妥", "舒降之", "美百乐镇",
            "血脂康", "脂必妥", "他汀", "贝特", "烟酸",
            # 补充
            "crestor", "lipitor", "zocor", "pravachol", "lescol", "mevacor",
            "vytorin", "zetia", "welchol", "questran", "lovaza", "trilipix",
            "pemafibrate", "pemfibrate",
        ],
        "category": "cardiovascular",
    },
    "甲状腺药物(Thyroid)": {
        "keywords": [
            "levothyroxine", "synthroid", "levoxyl", "euthyrox", "unithroid",
            "liothyronine", "cytomel", "armour thyroid", "nature-throid",
            "methimazole", "tapazole", "propylthiouracil", "ptu", "thyrolar",
            "左甲状腺素", "优甲乐", "雷替斯", "加衡", "泽宁",
            "甲巯咪唑", "赛治", "他巴唑", "丙硫氧嘧啶",
            "甲状腺片", "甲状腺素",
            # 补充
            "tirosint", "thyrogen", "levothroid", "np thyroid",
            "甲状旁腺",
        ],
        "category": "endocrine",
    },
    "降糖药(Diabetes)": {
        "keywords": [
            "metformin", "glipizide", "glyburide", "glimepiride", "gliclazide",
            "sitagliptin", "saxagliptin", "linagliptin", "alogliptin", "vildagliptin",
            "empagliflozin", "dapagliflozin", "canagliflozin", "ertugliflozin",
            "liraglutide", "semaglutide", "dulaglutide", "exenatide", "tirzepatide",
            "insulin", "novolog", "humalog", "lantus", "levemir", "tresiba",
            "pioglitazone", "rosiglitazone", "acarbose", "miglitol", "repaglinide",
            "nateglinide", "pramlintide", "diazoxide",
            "二甲双胍", "格列本脲", "格列美脲", "格列齐特", "格列吡嗪",
            "西格列汀", "沙格列汀", "利格列汀", "阿格列汀",
            "恩格列净", "达格列净", "卡格列净",
            "利拉鲁肽", "司美格鲁肽", "度拉糖肽",
            "胰岛素", "诺和灵", "优泌林", "来得时", "诺和锐",
            "阿卡波糖", "拜糖平", "格列酮", "瑞格列奈", "那格列奈",
            # 补充
            "januvia", "janumet", "glucophage", "actos", "avandia",
            "trulicity", "victoza", "ozempic", "rybelsus", "mounjaro",
            "farxiga", "invokana", "jardiance", "ste glat", "glyxambi",
            "synjardy", "trijardy", "diquil", "glucotrol", "micronase",
            "diabeta", "prandin", "starlix", "symlin",
            "艾塞那肽", "贝那鲁肽", "洛塞那肽",
        ],
        "category": "endocrine",
    },
    "降压药(Antihypertensive)": {
        "keywords": [
            "lisinopril", "enalapril", "ramipril", "captopril", "quinapril",
            "perindopril", "trandolapril", "benazepril", "fosinopril",
            "losartan", "valsartan", "irbesartan", "telmisartan", "candesartan",
            "olmesartan", "eprosartan", "azilsartan",
            "amlodipine", "nifedipine", "felodipine", "nicardipine",
            "diltiazem", "verapamil",
            "hydrochlorothiazide", "chlorthalidone", "furosemide", "bumetanide",
            "spironolactone", "eplerenone", "triamterene",
            "metoprolol", "atenolol", "propranolol", "carvedilol",
            "bisoprolol", "nebivolol", "labetalol",
            "doxazosin", "prazosin", "terazosin", "clonidine", "methyldopa",
            "hydralazine", "minoxidil",
            "普利", "沙坦", "地平", "洛尔", "噻嗪",
            "螺内酯", "呋塞米", "氢氯噻嗪",
            # 补充
            "cozaar", "hyzaar", "diovan", "exforge", "avalide",
            "norvasc", "cardizem", "lotrel", "tekturna", "mavik",
            "aceon", "accupril", "altace", "maxzide", "dyazide",
            "toprol", "lopressor", "tenormin", "coreg", "aldactone",
            "catapres", "loniten", "apresoline",
            "吲达帕胺", "寿比山", "硝苯地平", "拜新同", "络活喜",
            "代文", "科素亚", "安博维", "美卡素", "必洛斯",
            "倍他乐克", "博苏", "康可",
        ],
        "category": "cardiovascular",
    },
    "激素类(Steroids)": {
        "keywords": [
            "prednisone", "prednisolone", "dexamethasone", "betamethasone",
            "hydrocortisone", "cortisone", "methylprednisolone", "triamcinolone",
            "estrogen", "progesterone", "testosterone", "estradiol",
            "medroxyprogesterone", "norethindrone", "levonorgestrel",
            "conjugated estrogens", "premarin",
            "泼尼松", "泼尼松龙", "地塞米松", "甲泼尼龙",
            "氢化可的松", "倍他米松", "雌激素", "孕激素", "睾酮",
            "强的松", "可的松",
            # 补充
            "deltasone", "medrol", "solumedrol", "kenalog", "celestone",
            "depomedrol", "cortef", "florinef", "androgel", "axiron",
            "depo-provera", "ortho-evra", "nuvaring", "mirena",
            "雌二醇", "补佳乐", "倍美力", "黄体酮", "炔诺酮",
        ],
        "category": "endocrine",
    },
    "抗凝药(Anticoagulant)": {
        "keywords": [
            "warfarin", "coumadin", "heparin", "enoxaparin", "dalteparin",
            "rivaroxaban", "apixaban", "edoxaban", "dabigatran",
            "clopidogrel", "ticagrelor", "prasugrel", "aspirin",
            "dipyridamole", "ticlopidine", "cilostazol",
            "华法林", "阿司匹林", "氯吡格雷", "替格瑞洛",
            "利伐沙班", "阿哌沙班", "达比加群",
            "低分子肝素",
            # 补充
            "xarelto", "eliquis", "pradaxa", "savaysa", "plavix",
            "brilinta", "effient", "persantine", "aggrenox", "pletal",
            "lovenox", "fragmin", "innohep", "arixtra",
            "阿司匹林肠溶片", "拜阿司匹灵",
            "双嘧达莫", "潘生丁",
            "肝素钠", "依诺肝素",
        ],
        "category": "cardiovascular",
    },
}

# 疾病分类（用于排除既往病史）
DISEASE_EXCLUSION_CATEGORIES = {
    "心血管病(CVD)": {"vars": ["冠心病(CHD)", "心肌梗死(HeartAttack)", "心力衰竭(CHF)", "心绞痛(Angina)", "中风(Stroke)"]},
    "糖尿病(Diabetes)": {"vars": ["糖尿病诊断(DM-Dx)"]},
    "慢性肝病(LiverDz)": {"vars": ["慢性肝病(Liver)"]},
    "慢性肾病(CKD)": {"vars": ["慢性肾病(CKD)", "肾衰竭(KidneyFail)"]},
    "甲状腺病(ThyroidDz)": {"vars": ["甲状腺问题(Thyroid)", "甲状腺病(ThyroidDz)"]},
    "癌症(Cancer)": {"vars": ["癌症(Cancer)"]},
    "COPD/哮喘": {"vars": ["COPD", "哮喘(Asthma)", "慢性支气管炎(ChrBronch)", "肺气肿(Emphysema)"]},
    "关节炎/骨质疏松": {"vars": ["类风湿关节炎(RA)", "骨关节炎(Osteoarth)", "骨质疏松(Osteoporosis)", "痛风(Gout)"]},
    "痴呆(Dementia)": {"vars": ["痴呆/阿尔茨海默(Dementia)"]},
    "怀孕(Pregnant)": {"vars": ["怀孕状态(Pregnant)"]},
}

# ============================================================================
# 复合筛选条件系统（借鉴 数据规范化处理器.py 的 FilterCondition）
# ============================================================================
# 支持的操作符（与您已有工具完全一致）
FILTER_OPERATORS = ["≥", ">", "=", "≠", "<", "≤", "区间"]

# 在 profile JSON 中的格式：
# "compound_filters": {
#     "between_groups": "AND",    # 组间关系 AND/OR
#     "groups": [
#         {
#             "within_group": "AND",  # 组内关系 AND/OR
#             "label": "肝功能正常组",
#             "conditions": [
#                 {"var": "谷草转氨酶(AST,IU/L)", "op": "区间", "min": 3, "max": 50},
#                 {"var": "γ-谷氨酰转移酶(GGT,IU/L)", "op": "≤", "value": 40}
#             ]
#         },
#         {
#             "within_group": "AND",
#             "label": "血糖异常组",
#             "conditions": [
#                 {"var": "空腹血糖(Glu,mmol/L)", "op": "≥", "value": 7.0}
#             ]
#         }
#     ]
# }

def apply_compound_filter(df, filter_config):
    """
    应用复合筛选条件
    filter_config: {
        "between_groups": "AND"/"OR",
        "groups": [{
            "within_group": "AND"/"OR",
            "label": "...",
            "conditions": [{"var": ..., "op": ..., "min": ..., "max": ..., "value": ...}, ...]
        }, ...]
    }
    返回: 筛选后的DataFrame
    """
    if not filter_config or "groups" not in filter_config or not filter_config["groups"]:
        return df

    between_logic = filter_config.get("between_groups", "AND")
    group_results = []

    for group in filter_config["groups"]:
        within_logic = group.get("within_group", "AND")
        cond_results = []

        for cond in group.get("conditions", []):
            var_keyword = cond.get("var", "")
            op = cond.get("op", "=")
            value = cond.get("value")
            vmin = cond.get("min")
            vmax = cond.get("max")

            # 模糊查找列名
            found_col = None
            kw = var_keyword.upper().replace(" ", "")
            for c in df.columns:
                if kw in c.upper().replace(" ", ""):
                    found_col = c
                    break
            if not found_col:
                continue

            col_series = df[found_col]
            col_data = pd.to_numeric(col_series, errors="coerce")

            if op == "区间":
                if vmin is not None and vmax is not None:
                    cond_results.append(col_data.between(vmin, vmax))
            elif op == "≥":
                cond_results.append(col_data >= value)
            elif op == ">":
                cond_results.append(col_data > value)
            elif op == "=":
                cond_results.append(col_data == value)
            elif op == "≠":
                cond_results.append(col_data != value)
            elif op == "<":
                cond_results.append(col_data < value)
            elif op == "≤":
                cond_results.append(col_data <= value)

        if not cond_results:
            continue
        if within_logic == "AND":
            group_results.append(pd.concat(cond_results, axis=1).all(axis=1))
        else:
            group_results.append(pd.concat(cond_results, axis=1).any(axis=1))

    if not group_results:
        return df

    if between_logic == "AND":
        final_mask = pd.concat(group_results, axis=1).all(axis=1)
    else:
        final_mask = pd.concat(group_results, axis=1).any(axis=1)

    return df[final_mask].copy()


# ============================================================================
# 衍生变量计算规则（可扩展）
# ============================================================================
DERIVED_FORMULAS = {
    "TyG指数(TyG)": {
        "label": "TyG指数",
        "description": "胰岛素抵抗指数 = Ln(TG×Glu/2)",
        "requires": ["甘油三酯(TG,mg/dL)", "空腹血糖(Glu,mg/dL)"],
        "formula": lambda tg, glu: np.log(tg * glu / 2),
    },
    "TSH/FT4比值(TSH/FT4)": {
        "label": "TSH/FT4比值",
        "description": "甲状腺激素敏感性指数",
        "requires": ["促甲状腺激素(TSH,mIU/L)", "游离甲状腺素(FT4,ng/dL)"],
        "formula": lambda tsh, ft4: tsh / ft4,
    },
    "FT3/FT4比值(FT3/FT4)": {
        "label": "FT3/FT4比值",
        "description": "外周甲状腺激素转化效率",
        "requires": ["游离三碘甲腺原氨酸(FT3,pg/mL)", "游离甲状腺素(FT4,ng/dL)"],
        "formula": lambda ft3, ft4: ft3 / ft4,
    },
    "HOMA-IR": {
        "label": "HOMA-IR",
        "description": "胰岛素抵抗指数 = Glu×Ins/22.5 (Glu需mmol/L)",
        "requires": ["空腹血糖(Glu,mmol/L)", "胰岛素(Ins,uU/mL)"],
        "formula": lambda glu, ins: glu * ins / 22.5,
    },
    "动脉粥样硬化指数(AI)": {
        "label": "动脉粥样硬化指数(AI)",
        "description": "AI = (TC-HDL)/HDL",
        "requires": ["总胆固醇(TC,mmol/L)", "高密度脂蛋白(HDL,mmol/L)"],
        "formula": lambda tc, hdl: (tc - hdl) / hdl,
    },
    "非HDL-C(non-HDL)": {
        "label": "非HDL-C(mmol/L)",
        "description": "non-HDL = TC - HDL",
        "requires": ["总胆固醇(TC,mmol/L)", "高密度脂蛋白(HDL,mmol/L)"],
        "formula": lambda tc, hdl: tc - hdl,
    },
    "尿酸/肌酐比值(UA/Cr)": {
        "label": "UA/Cr比值",
        "description": "UA(mg/dL) / Cr(mg/dL)",
        "requires": ["尿酸(UA,mg/dL)", "肌酐(Cr,mg/dL)"],
        "formula": lambda ua, cr: ua / cr,
    },
    "BMI分级(BMI-Cat)": {
        "label": "BMI分级",
        "description": "中国标准 BMI 分类",
        "requires": ["BMI(kg/m²)"],
        "formula": None,  # 特殊处理
    },
    "eGFR(CKD-EPI)": {
        "label": "eGFR(mL/min/1.73m²)",
        "description": "CKD-EPI估算肾小球滤过率(需年龄性别肌酐)",
        "requires": ["肌酐(Cr,mg/dL)", "年龄(Age,岁)", "性别(Gender)"],
        "formula": None,  # 特殊处理
    },
    "血脂达标(LDL-C达标)": {
        "label": "血脂达标(LDL-C达标)",
        "description": "根据ASCVD风险分层判断LDL-C是否达标(<2.6mmol/L)",
        "requires": ["低密度脂蛋白(LDL,mmol/L)"],
        "formula": None,  # 特殊处理
    },
    "高血压分级(BP-Cat)": {
        "label": "高血压分级",
        "description": "中国高血压防治指南分类",
        "requires": ["收缩压-1(SBP1,mmHg)", "舒张压-1(DBP1,mmHg)"],
        "formula": None,  # 特殊处理
    },
}

# ============================================================================
# 临床正常值范围（用于标记异常值）
# ============================================================================
CLINICAL_RANGES = [
    {"var": "BMI(kg/m²)", "min": 12, "max": 60, "unit": "kg/m²"},
    {"var": "收缩压-1(SBP1,mmHg)", "min": 60, "max": 260},
    {"var": "舒张压-1(DBP1,mmHg)", "min": 30, "max": 160},
    {"var": "总胆固醇(TC,mmol/L)", "min": 1.0, "max": 15.0},
    {"var": "高密度脂蛋白(HDL,mmol/L)", "min": 0.2, "max": 4.0},
    {"var": "低密度脂蛋白(LDL,mmol/L)", "min": 0.2, "max": 10.0},
    {"var": "甘油三酯(TG,mmol/L)", "min": 0.1, "max": 15.0},
    {"var": "空腹血糖(Glu,mmol/L)", "min": 2.0, "max": 35.0},
    {"var": "促甲状腺激素(TSH,mIU/L)", "min": 0.01, "max": 100},
    {"var": "肌酐(Cr,μmol/L)", "min": 10, "max": 1500},
    {"var": "尿酸(UA,μmol/L)", "min": 50, "max": 900},
    {"var": "谷草转氨酶(AST,IU/L)", "min": 1, "max": 1000},
    {"var": "白细胞(WBC,×10⁹/L)", "min": 0.5, "max": 50},
    {"var": "血小板(PLT,×10⁹/L)", "min": 10, "max": 1000},
    {"var": "血红蛋白(Hb,g/L)", "min": 30, "max": 250},
]

# ============================================================================
# 研究方案系统（Study Profile）
# ============================================================================
PROFILES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "profiles")
os.makedirs(PROFILES_DIR, exist_ok=True)

DEFAULT_STUDY_PROFILE = {
    "study_name": "新研究方案",
    "description": "",
    "data": {
        "cycles": [],
        "groups": ["demo_core", "bmx", "bpx", "bpq", "smq", "alq", "diq",
                    "paq", "slq", "dpq", "huq", "lipid", "glu", "liver",
                    "kidney", "thy", "vit", "inf", "metal", "rxq", "diet"],
        "convert_units": True,
        "auto_aggregate": True,
        "include_mortality": False,
    },
    "inclusion": [],
    "exclusion": {
        "medication_classes": [],
        "range_filters": [],
        "disease_classes": [],
        "drug_audit": True,
    },
    "compound_filters": {},
    "derived_indices": [],
    "flag_outliers": True,
    "flag_missing_pct": 80,
}

STUDY_TEMPLATES = {
    "甲状腺功能研究": {
        "study_name": "甲状腺功能与代谢研究(J-Shaped TSH标准)",
        "description": "基于论文10项QC过滤器：TSH 0.1-10.0, FT3 0.5-50.0pmol/L, FT4 0.5-70.0pmol/L, TC 1.0-13.0, TG 0.2-11.3, HDL 0.2-3.0, LDL 0.3-8.0, ALT 3-200, Cr 10-300μmol/L, 年龄≥18, BMI 16-40",
        "data": {"cycles": ["E","F","G"], "groups": ["demo_core","demo_detail","bmx","bpx","bpq",
                 "smq","alq","diq","lipid","glu","liver","kidney","thy","inf",
                 "mcq","rxq"],
                 "convert_units": True, "auto_aggregate": True},
        "inclusion": [
            {"var": "年龄(Age,岁)", "min": 18, "reason": "年龄≥18岁"},
            {"var": "BMI(kg/m²)", "min": 16, "max": 40, "reason": "BMI 16-40 kg/m²"},
        ],
        "exclusion": {
            "medication_classes": ["甲状腺药物(Thyroid)"],
            "range_filters": [
                {"var": "促甲状腺激素(TSH,mIU/L)", "min": 0.1, "max": 10.0, "reason": "TSH 0.1-10.0 mIU/L"},
                {"var": "游离三碘甲腺原氨酸(FT3,pmol/L)", "min": 0.5, "max": 50.0, "reason": "FT3 0.5-50.0 pmol/L"},
                {"var": "游离甲状腺素(FT4,pmol/L)", "min": 0.5, "max": 70.0, "reason": "FT4 0.5-70.0 pmol/L"},
                {"var": "总胆固醇(TC,mmol/L)", "min": 1.0, "max": 13.0, "reason": "TC 1.0-13.0 mmol/L"},
                {"var": "甘油三酯(TG,mmol/L)", "min": 0.2, "max": 11.3, "reason": "TG 0.2-11.3 mmol/L"},
                {"var": "高密度脂蛋白(HDL,mmol/L)", "min": 0.2, "max": 3.0, "reason": "HDL 0.2-3.0 mmol/L"},
                {"var": "低密度脂蛋白(LDL,mmol/L)", "min": 0.3, "max": 8.0, "reason": "LDL 0.3-8.0 mmol/L"},
                {"var": "谷丙转氨酶(ALT,IU/L)", "min": 3, "max": 200, "reason": "ALT 3-200 U/L"},
                {"var": "肌酐(Cr,μmol/L)", "min": 10, "max": 300, "reason": "Cr 10-300 μmol/L"},
            ],
            "disease_classes": [],
            "drug_audit": True,
        },
        "derived_indices": ["TSH/FT4比值(TSH/FT4)", "FT3/FT4比值(FT3/FT4)", "TyG指数(TyG)",
                            "动脉粥样硬化指数(AI)", "非HDL-C(non-HDL)", "血脂达标(LDL-C达标)",
                            "高血压分级(BP-Cat)"],
        "flag_outliers": True, "flag_missing_pct": 80,
    },
    "代谢综合征研究": {
        "study_name": "代谢综合征影响因素分析",
        "description": "纳入全部可用周期，排除降糖药、降脂药使用者及心血管病史",
        "data": {"cycles": ["H","I","J"], "groups": ["demo_core","bmx","bpx","smq",
                 "alq","diq","lipid","glu","ins","liver","kidney","inf","mcq","rxq"],
                 "convert_units": True, "auto_aggregate": True},
        "inclusion": [{"var": "年龄(Age,岁)", "min": 20, "max": 80}],
        "exclusion": {
            "medication_classes": ["降糖药(Diabetes)", "降脂药(Lipid-lowering)"],
            "range_filters": [],
            "disease_classes": ["心血管病(CVD)"],
            "drug_audit": True,
        },
        "derived_indices": ["TyG指数(TyG)", "HOMA-IR", "动脉粥样硬化指数(AI)", "非HDL-C(non-HDL)"],
        "flag_outliers": True, "flag_missing_pct": 70,
    },
    "全因死亡率研究": {
        "study_name": "全因死亡率影响因素",
        "description": "需要选择任意周期并勾选右侧'包含死亡率数据'",
        "data": {"cycles": ["H","I","J"], "groups": ["demo_core","bmx","bpx","smq",
                 "alq","diq","lipid","glu","liver","kidney","inf"],
                 "convert_units": True, "auto_aggregate": True, "include_mortality": True},
        "inclusion": [{"var": "年龄(Age,岁)", "min": 18}],
        "exclusion": {"medication_classes": [], "range_filters": []},
        "derived_indices": [],
        "flag_outliers": True, "flag_missing_pct": 80,
    },
    "Wang2023复刻(甲状腺+睡眠)": {
        "study_name": "Wang 2023 精确复刻 — PLoS One 2023",
        "description": "精确复刻Wang 2023 (PLoS One)筛选：年龄≥18 + 排除甲状腺病/药/怀孕 + 有睡眠数据。不加TSH范围限制",
        "data": {"cycles": ["E","F","G"], "groups": ["demo_core","demo_detail","bmx","bpx","bpq",
                 "smq","alq","diq","lipid","glu","liver","kidney","thy","inf",
                 "mcq","rxq","slq"],
                 "convert_units": True, "auto_aggregate": True},
        "inclusion": [{"var": "年龄(Age,岁)", "min": 18}],
        "exclusion": {
            "medication_classes": ["甲状腺药物(Thyroid)"],
            "disease_classes": ["甲状腺病(ThyroidDz)", "怀孕(Pregnant)"],
            "range_filters": [],   # Wang不限制TSH/FT3/FT4范围
            "drug_audit": True,
        },
        "derived_indices": ["TSH/FT4比值(TSH/FT4)", "FT3/FT4比值(FT3/FT4)", "TyG指数(TyG)",
                            "动脉粥样硬化指数(AI)", "非HDL-C(non-HDL)", "血脂达标(LDL-C达标)",
                            "高血压分级(BP-Cat)"],
        "flag_outliers": True, "flag_missing_pct": 80,
    },
}

class StudyProfile:
    """研究方案配置：加载/保存/应用研究方案"""

    def __init__(self, profile_path=None):
        self.path = profile_path
        self.config = dict(DEFAULT_STUDY_PROFILE)

    def load(self, path=None):
        fp = path or self.path
        if not fp or not os.path.exists(fp):
            return False
        try:
            with open(fp, "r", encoding="utf-8") as f:
                loaded = json.load(f)
            # 合并保留默认值
            merged = dict(DEFAULT_STUDY_PROFILE)
            merged.update(loaded)
            if "data" in loaded:
                merged["data"].update(loaded["data"])
            if "exclusion" in loaded:
                merged["exclusion"].update(loaded["exclusion"])
            self.config = merged
            self.path = fp
            return True
        except Exception as e:
            print(f"加载方案失败: {e}")
            return False

    def save(self, path=None):
        fp = path or self.path
        if not fp:
            fp = os.path.join(PROFILES_DIR, f"{self.config['study_name']}.json")
        try:
            with open(fp, "w", encoding="utf-8") as f:
                json.dump(self.config, f, ensure_ascii=False, indent=2)
            self.path = fp
            return True
        except Exception as e:
            print(f"保存方案失败: {e}")
            return False

    @staticmethod
    def list_profiles():
        """列出所有可用方案"""
        profiles = []
        for fname in os.listdir(PROFILES_DIR):
            if fname.endswith(".json"):
                fpath = os.path.join(PROFILES_DIR, fname)
                try:
                    with open(fpath, "r", encoding="utf-8") as f:
                        cfg = json.load(f)
                    profiles.append({
                        "path": fpath,
                        "name": cfg.get("study_name", fname),
                        "desc": cfg.get("description", ""),
                    })
                except:
                    profiles.append({"path": fpath, "name": fname, "desc": ""})
        return profiles

    def get_selected_cycles(self):
        """获取方案配置的周期后缀"""
        return list(self.config["data"].get("cycles", []))

    def get_selected_groups(self):
        """获取方案配置的变量组"""
        return list(self.config["data"].get("groups", []))

    def get_params(self):
        """获取引擎运行参数"""
        d = self.config["data"]
        return {
            "selected_cycle_suffixes": d.get("cycles", []),
            "selected_group_keys": d.get("groups", []),
            "include_mortality": d.get("include_mortality", False),
            "convert_units": d.get("convert_units", True),
            "auto_aggregate": d.get("auto_aggregate", True),
        }


class LogReporter:
    """筛选流水账记录器：记录每步筛选的筛前/筛后人数，最终输出TXT报告"""

    def __init__(self):
        self.records = []  # [(step_name, before, after, reason), ...]
        self.start_time = datetime.now()

    def record(self, step_name, before, after, reason=""):
        """记录一条筛选步骤"""
        self.records.append((step_name, before, after, reason))

    def write_txt(self, output_path, final_rows, final_cols, cycle_dist=None):
        """生成TXT筛选报告，与输出文件放在同一目录"""
        if not self.records:
            return

        # 生成报告路径: 与输出文件同目录，文件名+_筛选报告.txt
        base_dir = os.path.dirname(output_path) if output_path else os.getcwd()
        base_name = os.path.splitext(os.path.basename(output_path))[0] if output_path else "NHANES"
        report_path = os.path.join(base_dir, f"{base_name}_筛选报告.txt")

        elapsed = (datetime.now() - self.start_time).total_seconds()

        lines = []
        lines.append("=" * 80)
        lines.append(f"                NHANES 数据筛选清洗报告")
        lines.append(f"                生成时间: {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append(f"                耗时: {elapsed:.1f} 秒")
        lines.append("=" * 80)
        lines.append("")

        # 一、纳入标准
        inclusion = [r for r in self.records if "纳入" in r[0]]
        exclusion = [r for r in self.records if "纳入" not in r[0]]

        lines.append("一、纳入标准")
        lines.append("─" * 80)
        if inclusion:
            for i, (name, before, after, reason) in enumerate(inclusion, 1):
                removed = before - after
                lines.append(f"  {i}. {name}")
                if reason:
                    lines.append(f"     条件: {reason}")
                lines.append(f"     筛前: {before:,} → 筛后: {after:,} → 排除: {removed:,} 人")
                lines.append("")
        else:
            lines.append("  (无纳入条件)")
            lines.append("")

        lines.append("二、排除标准")
        lines.append("─" * 80)
        if exclusion:
            for i, (name, before, after, reason) in enumerate(exclusion, 1):
                removed = before - after
                lines.append(f"  {i}. {name}")
                if reason:
                    lines.append(f"     条件: {reason}")
                lines.append(f"     筛前: {before:,} → 筛后: {after:,} → 排除: {removed:,} 人")
                lines.append("")
        else:
            lines.append("  (无排除条件)")
            lines.append("")

        # 最终汇总
        lines.append("三、最终数据汇总")
        lines.append("─" * 80)
        lines.append(f"  总记录数: {final_rows:,} 人")
        lines.append(f"  总变量数: {final_cols} 列")
        if cycle_dist:
            lines.append(f"  周期分布: {cycle_dist}")
        lines.append("")
        
        # 四、偏态变量警告（新增）
        if hasattr(self, 'qc_report') and self.qc_report.get("偏态变量"):
            skewed_vars = self.qc_report["偏态变量"]
            if skewed_vars:
                lines.append("四、⚠️ 偏态变量警告（建议报告中位数而非均值）")
                lines.append("─" * 80)
                lines.append("  【偏态定义】偏度>1为右偏态，偏度<-1为左偏态")
                lines.append("  【报告建议】偏态变量应报告：中位数(IQR)，而非均值±标准差")
                lines.append("  【统计建议】参数检验前需对数转换，否则建议用非参数检验")
                lines.append("")
                for var, info in skewed_vars.items():
                    lines.append(f"  ▸ {var}")
                    lines.append(f"    偏度: {info['偏度']} ({info['类型']})")
                    lines.append(f"    中位数(IQR): {info['中位数']} ({info['IQR']})")
                    lines.append(f"    建议: {info['建议']}")
                    lines.append("")
        
        lines.append("=" * 80)
        lines.append("  报告结束")
        lines.append("=" * 80)

        content = "\n".join(lines)
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(content)

        return report_path

    def get_summary_lines(self):
        """返回每步筛选的摘要文本（用于嵌入CSV文件头）"""
        lines = []
        for name, before, after, reason in self.records:
            removed = before - after
            if reason:
                lines.append(f"筛选: {name} | 条件: {reason} | {before:,}→{after:,} | 排除:{removed:,}")
            else:
                lines.append(f"筛选: {name} | {before:,}→{after:,} | 排除:{removed:,}")
        return lines


class CleaningPipeline:
    """清洗管道：纳入排除 + 药物排除 + 衍生变量 + 异常值标记"""

    def __init__(self, df, profile: StudyProfile, rxq_data=None, reporter=None):
        self.df = df.copy()
        self.profile = profile
        self.rxq_data = rxq_data  # 可选：RXQ_RX 表用于药物排除
        self.reporter = reporter  # LogReporter 实例
        self.logs = []
        self.removed_count = 0
        self.columns = set(df.columns)

    def log(self, msg):
        self.logs.append(msg)
        print(f"  [清洗] {msg}")

    def _find_col(self, keyword):
        """模糊查找列名，支持精确→子串→分词三级回退"""
        kw = keyword.upper().replace(" ", "")
        # 1级：精确子串匹配
        for c in self.columns:
            if kw in c.upper().replace(" ", ""):
                return c
        # 2级：更宽松的匹配（取逗号前第一部分）
        kw_parts = keyword.split(",")[0].strip().upper()
        for c in self.columns:
            if kw_parts in c.upper():
                return c
        # 3级：分词匹配——英文变量名每个词分别匹配
        kw_words = [w for w in kw.replace("(", " ").replace(")", " ").replace("/", " ").split() if len(w) > 2]
        if kw_words:
            for c in self.columns:
                c_up = c.upper()
                if all(any(w in c_up_part for c_up_part in [c_up]) or
                       any(w in part for part in c_up.split()) for w in kw_words[:3]):
                    # 至少匹配2个词才算
                    match_count = sum(1 for w in kw_words if w in c_up)
                    if match_count >= 2:
                        return c
        return None

    def run(self):
        """执行完整清洗管道"""
        before = len(self.df)
        cfg = self.profile.config

        # 1. 纳入标准
        for criterion in cfg.get("inclusion", []):
            self._apply_inclusion(criterion)

        # 2. 排除标准 - 范围过滤
        for filt in cfg["exclusion"].get("range_filters", []):
            self._apply_exclusion(filt)

        # 3. 排除标准 - 药物排除（含审核报告）
        med_classes = cfg["exclusion"].get("medication_classes", [])
        if med_classes:
            self._apply_medication_exclusion(med_classes)
            # 生成药物审核报告
            drug_audit = cfg.get("drug_audit", False) or cfg.get("flag_outliers", True)
            if drug_audit:
                self._generate_drug_audit_report()

        # 3b. 排除标准 - 疾病排除
        disease_classes = cfg["exclusion"].get("disease_classes", [])
        if disease_classes:
            self._apply_disease_exclusion(disease_classes)

        # 3c. 复合筛选条件（支持 AND/OR 多组条件）
        compound_filters = cfg.get("compound_filters")
        if compound_filters and compound_filters.get("groups"):
            before_f = len(self.df)
            self.df = apply_compound_filter(self.df, compound_filters)
            after_f = len(self.df)
            removed_f = before_f - after_f
            if removed_f > 0:
                n_groups = len(compound_filters["groups"])
                self.log(f"✅ 复合条件筛选: 排除 {removed_f} 人 ({n_groups} 组条件)")
                if self.reporter:
                    self.reporter.record(f"复合条件筛选({n_groups}组)", before_f, after_f,
                                         f"{n_groups} 组复合条件")

        # 4. 衍生变量
        for idx_name in cfg.get("derived_indices", []):
            self._apply_derived(idx_name)

        # 5. 异常值标记
        if cfg.get("flag_outliers", True):
            self._flag_outliers()

        # 6. 高缺失率变量标记
        max_miss = cfg.get("flag_missing_pct", 80)
        self._flag_high_missing(max_miss)

        after = len(self.df)
        if before > after:
            self.log(f"纳入排除: {before} → {after} (排除 {before-after} 人)")
        self.removed_count = before - after
        return self.df

    def _apply_inclusion(self, criterion):
        """应用纳入标准"""
        col = self._find_col(criterion.get("var", ""))
        if not col:
            self.log(f"⏭️ 纳入条件未找到列: {criterion['var']}")
            return
        vmin = criterion.get("min")
        vmax = criterion.get("max")
        not_null = criterion.get("not_null", False)
        before = len(self.df)
        mask = pd.Series(True, index=self.df.index)
        if vmin is not None:
            mask &= pd.to_numeric(self.df[col], errors="coerce") >= vmin
        if vmax is not None:
            mask &= pd.to_numeric(self.df[col], errors="coerce") <= vmax
        if not_null:
            mask &= self.df[col].notna()
        self.df = self.df[mask].copy()
        removed = before - len(self.df)
        if removed > 0:
            if not_null:
                reason = criterion.get("reason", f"{criterion['var']} 非空")
            else:
                reason = criterion.get("reason", f"{criterion['var']} {vmin}~{vmax}")
            self.log(f"✅ 纳入筛选: {col} → 排除 {removed} 人 ({reason})")
            if self.reporter:
                self.reporter.record(f"纳入筛选: {col}", before, len(self.df), reason)

    def _apply_exclusion(self, filt):
        """应用范围排除（自动匹配单位）"""
        col = self._find_col(filt.get("var", ""))
        if not col:
            self.log(f"⚠️ 范围排除: 列 '{filt.get('var', '')}' 在当前周期数据中不存在，跳过")
            return
        vmin = filt.get("min")
        vmax = filt.get("max")
        
        # ---- 单位自动匹配 ----
        # 提取指定列名中的单位（如"总胆固醇(TC,mmol/L)"→"mmol/L"）
        spec_unit = self._extract_unit(filt.get("var", ""))
        # 提取实际列名中的单位（如"总胆固醇(TC,mg/dL)"→"mg/dL"）
        col_unit = self._extract_unit(col)
        if spec_unit and col_unit and spec_unit != col_unit:
            # 转换阈值到列的单位
            converted = self._convert_threshold(vmin, vmax, spec_unit, col_unit, filt.get("var", ""))
            if converted:
                vmin, vmax = converted
                self.log(f"  ↪ 自动单位转换: {spec_unit} -> {col_unit} (阈值已调整)")
        
        before = len(self.df)
        mask = pd.Series(True, index=self.df.index)
        if vmin is not None:
            mask &= pd.to_numeric(self.df[col], errors="coerce") >= vmin
        if vmax is not None:
            mask &= pd.to_numeric(self.df[col], errors="coerce") <= vmax
        self.df = self.df[mask].copy()
        removed = before - len(self.df)
        if removed > 0:
            reason = filt.get("reason", f"{col} 范围排除")
            self.log(f"✅ 范围排除: {col} → 排除 {removed} 人 ({reason})")
            if self.reporter:
                self.reporter.record(f"范围排除: {col}", before, len(self.df), reason)

    @staticmethod
    def _extract_unit(col_name):
        """从列名中提取单位，如 '总胆固醇(TC,mg/dL)' → 'mg/dL'"""
        import re
        m = re.search(r',([a-zA-Z/μ²³]+)\)', col_name)
        if m:
            return m.group(1).lower().replace('μ', 'u').replace('²', '2').replace('³', '3')
        return ""

    UNIT_MAP = {
        "总胆固醇": ("mmol/l", "mg/dl", 38.67),
        "高密度脂蛋白": ("mmol/l", "mg/dl", 38.67),
        "低密度脂蛋白": ("mmol/l", "mg/dl", 38.67),
        "甘油三酯": ("mmol/l", "mg/dl", 88.57),
        "空腹血糖": ("mmol/l", "mg/dl", 18.01),
        "肌酐": ("umol/l", "mg/dl", 88.4),
        "游离甲状腺素": ("pmol/l", "ng/dl", 12.87),
        "游离三碘甲腺原氨酸": ("pmol/l", "pg/ml", 1.536),
    }

    @classmethod
    def _convert_threshold(cls, vmin, vmax, from_unit, to_unit, var_spec=""):
        """如果列的单位与设定的单位不同，转换阈值"""
        for var_name, (src_u, dst_u, factor) in cls.UNIT_MAP.items():
            if var_name not in var_spec:
                continue
            if (from_unit == src_u and to_unit == dst_u) or \
               (from_unit == dst_u and to_unit == src_u):
                f = factor if from_unit == src_u else 1.0 / factor
                return (
                    vmin * f if vmin is not None else None,
                    vmax * f if vmax is not None else None,
                )
        return None

    def _apply_medication_exclusion(self, med_classes):
        """排除服用特定药物的参与者"""
        if self.rxq_data is None:
            err_msg = "❌ 需要药物排除但未加载 RXQ_RX 数据！请在左侧变量组中勾选「💊 处方药」后重试"
            self.log(err_msg)
            raise ValueError(err_msg)
        if "RXDDRUG" not in self.rxq_data.columns:
            err_msg = "❌ RXQ_RX 表中无 RXDDRUG 列，无法进行药物排除"
            self.log(err_msg)
            raise ValueError(err_msg)

        excluded_seqns = set()
        drug_names = self.rxq_data["RXDDRUG"].dropna().astype(str).str.lower()
        for cls_name in med_classes:
            if cls_name not in DRUG_CLASSIFICATION:
                continue
            keywords = DRUG_CLASSIFICATION[cls_name]["keywords"]
            mask = pd.Series(False, index=drug_names.index)
            for kw in keywords:
                mask |= drug_names.str.contains(kw.lower(), na=False)
            matched = self.rxq_data[mask]["SEQN"].unique()
            excluded_seqns.update(matched)
            self.log(f"✅ {cls_name}: 匹配 {len(matched)} 人")

        seqn_col = self._find_col("序号(SEQN)")
        if seqn_col and excluded_seqns:
            before = len(self.df)
            seqn_vals = set(pd.to_numeric(self.df[seqn_col], errors="coerce"))
            to_remove = seqn_vals & excluded_seqns
            self.df = self.df[~pd.to_numeric(self.df[seqn_col], errors="coerce").isin(to_remove)].copy()
            removed = before - len(self.df)
            self.log(f"✅ 药物排除: 共排除 {removed} 人 (涉及 {len(med_classes)} 类药物)")
            if self.reporter:
                self.reporter.record(f"药物排除: {', '.join(med_classes)}", before, len(self.df),
                                     f"涉及 {len(med_classes)} 类药物")

    def _apply_derived(self, idx_name):
        """计算衍生变量"""
        if idx_name not in DERIVED_FORMULAS:
            return
        formula = DERIVED_FORMULAS[idx_name]
        req_cols = []
        for req in formula["requires"]:
            found = self._find_col(req)
            if found:
                req_cols.append(found)
            else:
                self.log(f"⏭️ 衍生变量 {idx_name}: 缺少 {req}")
                return

        label = formula["label"]
        if label in self.df.columns:
            return  # 已存在

        try:
            if idx_name == "BMI分级(BMI-Cat)":
                col = req_cols[0]
                bmi = pd.to_numeric(self.df[col], errors="coerce")
                cat = pd.cut(bmi, bins=[0, 18.5, 24, 28, 100],
                             labels=["偏瘦", "正常", "超重", "肥胖"])
                self.df[label] = cat
                self.log(f"✅ 衍生: {label}")
            elif idx_name == "eGFR(CKD-EPI)":
                cr_col, age_col, gender_col = req_cols
                cr = pd.to_numeric(self.df[cr_col], errors="coerce")
                age = pd.to_numeric(self.df[age_col], errors="coerce")
                gender = self.df[gender_col].astype(str)
                # 简化MDRD公式
                is_female = gender.str.contains("女|Female", na=False)
                egfr = 175 * (cr/88.4)**(-1.154) * age**(-0.203)
                egfr = egfr * 0.742 if is_female.any() else egfr
                self.df[label] = egfr.round(0)
                self.log(f"✅ 衍生: {label}")
            elif idx_name == "血脂达标(LDL-C达标)":
                col = req_cols[0]
                ldl = pd.to_numeric(self.df[col], errors="coerce")
                self.df[label] = ldl.apply(lambda x: "达标" if pd.notna(x) and x < 2.6 else
                                           ("未达标" if pd.notna(x) and x >= 2.6 else np.nan))
                self.log(f"✅ 衍生: {label}")
            elif idx_name == "高血压分级(BP-Cat)":
                sbp_col, dbp_col = req_cols
                sbp = pd.to_numeric(self.df[sbp_col], errors="coerce")
                dbp = pd.to_numeric(self.df[dbp_col], errors="coerce")
                cat = pd.Series(np.nan, index=self.df.index)
                # 中国高血压防治指南分类
                cat.loc[(sbp < 120) & (dbp < 80)] = "正常"
                cat.loc[((sbp >= 120) & (sbp < 140)) | ((dbp >= 80) & (dbp < 90))] = "正常高值"
                cat.loc[((sbp >= 140) & (sbp < 160)) | (dbp >= 90)] = "高血压1级"
                cat.loc[(sbp >= 160) | (dbp >= 100)] = "高血压2级"
                cat.loc[(sbp >= 180) | (dbp >= 110)] = "高血压3级"
                self.df[label] = cat
                self.log(f"✅ 衍生: {label}")
            elif formula["formula"] is not None:
                vals = [pd.to_numeric(self.df[c], errors="coerce") for c in req_cols]
                result = formula["formula"](*vals)
                self.df[label] = result.round(3)
                self.log(f"✅ 衍生: {label}")
        except Exception as e:
            self.log(f"⚠️ 衍生 {idx_name} 失败: {e}")

    def _flag_outliers(self):
        """标记超出临床合理范围的值"""
        flag_cols = []
        for rule in CLINICAL_RANGES:
            col = self._find_col(rule["var"])
            if not col:
                continue
            vmin = rule.get("min")
            vmax = rule.get("max")
            mask = pd.Series(False, index=self.df.index)
            if vmin is not None:
                mask |= pd.to_numeric(self.df[col], errors="coerce") < vmin
            if vmax is not None:
                mask |= pd.to_numeric(self.df[col], errors="coerce") > vmax
            n_flag = mask.sum()
            if n_flag > 0:
                flag_name = f"⚠️{col}_异常"
                if flag_name not in self.df.columns:
                    self.df[flag_name] = mask.map({True: "异常", False: ""})
                    flag_cols.append(flag_name)
        if flag_cols:
            self.log(f"✅ 异常值标记: {len(flag_cols)} 个变量有异常标记")

    def _flag_high_missing(self, max_pct):
        """标记高缺失率变量"""
        flagged = []
        for col in self.df.columns:
            if col.startswith("⚠️") or col in ("序号(SEQN)", "调查周期"):
                continue
            pct = self.df[col].isna().mean() * 100
            if pct > max_pct:
                flagged.append(f"{col}({pct:.0f}%)")
        if flagged:
            self.log(f"⚠️ 高缺失率(>{max_pct}%): {', '.join(flagged[:5])}{'...' if len(flagged)>5 else ''}")

    def _apply_disease_exclusion(self, disease_classes):
        """根据既往病史排除参与者"""
        for cls_name in disease_classes:
            if cls_name not in DISEASE_EXCLUSION_CATEGORIES:
                self.log(f"⏭️ 未知疾病分类: {cls_name}")
                continue
            var_names = DISEASE_EXCLUSION_CATEGORIES[cls_name]["vars"]
            matched_cols = []
            for vn in var_names:
                col = self._find_col(vn)
                if col:
                    matched_cols.append(col)
            if not matched_cols:
                # 检查这些变量对应的原始变量名是否在当前数据源中存在
                raw_var_hints = []
                for vn in var_names:
                    # 尝试反向查找原始变量名
                    for g in VARIABLE_GROUPS:
                        for v in g["variables"]:
                            if v["var_label"] == vn:
                                raw_var_hints.append(v["var_name"])
                                break
                hint_str = f" (原始变量: {', '.join(raw_var_hints)})" if raw_var_hints else ""
                self.log(f"⏭️ {cls_name}: 当前周期数据表中无{var_names}这些列{hint_str}，跳过疾病排除")
                continue
            before = len(self.df)
            mask = pd.Series(False, index=self.df.index)
            for col in matched_cols:
                yes_vals = ["是", "1", 1, "Yes"]
                mask |= self.df[col].isin(yes_vals)
            self.df = self.df[~mask].copy()
            removed = before - len(self.df)
            if removed:
                self.log(f"✅ 疾病排除 {cls_name}: 排除 {removed} 人 (变量: {', '.join(matched_cols)})")
                if self.reporter:
                    self.reporter.record(f"疾病排除: {cls_name}", before, len(self.df),
                                         f"变量: {', '.join(matched_cols)}")

    def _generate_drug_audit_report(self):
        """生成药物排除审核报告（显示匹配了哪些药、漏了哪些药）"""
        if self.rxq_data is None or "RXDDRUG" not in self.rxq_data.columns:
            return
        all_drugs = self.rxq_data["RXDDRUG"].dropna().astype(str).str.lower().unique()
        med_classes = self.profile.config["exclusion"].get("medication_classes", [])
        matched_drugs = set()
        unmatched_drugs = []
        for drug in all_drugs:
            matched = False
            for cls_name in med_classes:
                if cls_name not in DRUG_CLASSIFICATION:
                    continue
                for kw in DRUG_CLASSIFICATION[cls_name]["keywords"]:
                    if kw.lower() in drug:
                        matched = True
                        matched_drugs.add(f"{drug[:40]}... → {cls_name}" if len(drug) > 40 else f"{drug} → {cls_name}")
                        break
                if matched:
                    break
            if not matched and len(drug) > 3:
                unmatched_drugs.append(drug)
        self.log(f"📋 药物审核: 匹配 {len(matched_drugs)} 种药名，未匹配 {len(unmatched_drugs)} 种")
        if matched_drugs:
            for md in sorted(list(matched_drugs))[:10]:
                self.log(f"   ✅ {md}")
            if len(matched_drugs) > 10:
                self.log(f"   ... 还有 {len(matched_drugs)-10} 种已匹配药名")
        if unmatched_drugs:
            self.log(f"   ⚠️ 以下药物未被任何分类覆盖（可能有遗漏）:")
            for ud in sorted(unmatched_drugs)[:15]:
                self.log(f"      ? {ud[:50]}")
            if len(unmatched_drugs) > 15:
                self.log(f"      ... 还有 {len(unmatched_drugs)-15} 种")
            self.log(f"   💡 如需补充药物，编辑 DRUG_CLASSIFICATION 字典或联系作者")


class NhanesGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("NHANES 数据下载与处理工具 v2.0")
        self.root.geometry("1200x850")
        self.root.minsize(1000, 700)

        self.engine = NhanesEngine()
        self.engine.log_callback = self._log

        self.cycle_vars = {}
        self.group_vars = {}
        self.last_result = None
        self.last_db_result = None
        self.db_builder = None
        self.db_import_btn = None
        self.custom_templates = []
        self._loaded_template_name = None  # 当前加载的模板名, 用于range_filters

        self._build_ui()
        self._check_db_status_on_startup()  # 检查本地数据库是否存在
        self._log("NHANES 数据工具 v2.0 已启动")
        self._log("")
        self._log("💡 使用说明:")
        self._log("  1. 左侧选择调查周期(年份)")
        self._log("  2. 中间选择SCI指标组，或在底部输入自定义变量")
        self._log("  3. 右上角可选「死亡率数据」")
        self._log("  4. 点击「开始处理」")
        self._log("")

    # noinspection PyAttributeOutsideInit
    def _build_ui(self):
        main = ttk.Frame(self.root, padding="8")
        main.pack(fill="both", expand=True)

        # 标题
        title_f = ttk.Frame(main)
        title_f.pack(fill="x")
        Label(title_f, text="🏥 NHANES 数据下载与处理工具 v2.0",
              font=("微软雅黑", 16, "bold"), fg="#2c3e50").pack(side="left")
        Label(title_f, text="美国国家健康与营养调查  |  CDC数据直连",
              font=("微软雅黑", 9), fg="#7f8c8d").pack(side="left", padx=15)

        body = ttk.PanedWindow(main, orient="horizontal")
        body.pack(fill="both", expand=True, pady=5)

        # ==== 左: 周期 ====
        left = ttk.LabelFrame(body, text="📅 调查周期", padding="5")
        body.add(left, weight=1)
        self._build_cycle_panel(left)

        # ==== 中: 变量（5x5网格，直接可见）====
        mid = ttk.LabelFrame(body, text="📊 数据指标组（25组 5×5排列）", padding="5")
        body.add(mid, weight=5)
        self._build_variable_panel(mid)

        # ==== 右: 控制 ====
        right = ttk.LabelFrame(body, text="⚙️ 操作", padding="5")
        body.add(right, weight=1)
        self._build_control_panel(right)

        # ==== 日志 ====
        log_f = ttk.LabelFrame(main, text="📝 运行日志", padding="3")
        log_f.pack(fill="both", pady=(5, 0))
        self._build_log_panel(log_f)

        # 进度
        self.progress = ttk.Progressbar(main, mode="determinate")
        self.progress.pack(fill="x", pady=(3, 0))

    def _build_cycle_panel(self, parent):
        btn_f = ttk.Frame(parent)
        btn_f.pack(fill="x", pady=2)
        ttk.Button(btn_f, text="全选", width=5,
                   command=lambda: [v.set(True) for v in self.cycle_vars.values()]).pack(side="left", padx=1)
        ttk.Button(btn_f, text="清空", width=5,
                   command=lambda: [v.set(False) for v in self.cycle_vars.values()]).pack(side="left", padx=1)

        for c_n, c_y, c_s in NHANES_CYCLES:
            v = BooleanVar(value=False)
            self.cycle_vars[c_s] = v
            ttk.Checkbutton(parent, text=f"{c_n}", variable=v).pack(anchor="w", padx=5, pady=1)

    def _build_variable_panel(self, parent):
        """直接显示5x5网格，下方加笔记本存放自定义变量/浏览器/方案"""
        # 主要区域：5x5网格（占大部分空间）
        grid_frame = ttk.Frame(parent)
        grid_frame.pack(fill="both", expand=True)
        self._build_preset_grid(grid_frame)

        # 下方：小号笔记本存放其他功能
        mini_book = ttk.Notebook(parent)
        mini_book.pack(fill="x", pady=(3, 0))

        # Tab: 自定义变量
        custom_frame = ttk.Frame(mini_book)
        mini_book.add(custom_frame, text="自定义变量")
        self._build_custom_tab(custom_frame)

        # Tab: 变量浏览器
        browser_frame = ttk.Frame(mini_book)
        mini_book.add(browser_frame, text="变量浏览器")
        self._build_browser_tab(browser_frame)

        # Tab: 研究方案
        profile_frame = ttk.Frame(mini_book)
        mini_book.add(profile_frame, text="研究方案")
        self._build_profile_tab(profile_frame)

    def _build_preset_grid(self, parent):
        """25个变量组排成5×5网格，无滚动条"""
        # 控制条
        ctrl = ttk.Frame(parent)
        ctrl.pack(fill="x", pady=(0, 3))
        ttk.Button(ctrl, text="全选", width=8,
                   command=lambda: [v.set(True) for v in self.group_vars.values()]).pack(side="left", padx=1)
        ttk.Button(ctrl, text="清空", width=8,
                   command=lambda: [v.set(False) for v in self.group_vars.values()]).pack(side="left", padx=1)
        ttk.Button(ctrl, text="推荐SCI指标", width=14,
                   command=self._select_sci).pack(side="left", padx=5)
        Label(ctrl, text="每组括号内为变量数", font=("微软雅黑", 7), fg="#999").pack(side="left", padx=5)

        # 5×5 网格
        grid = ttk.Frame(parent)
        grid.pack(fill="both", expand=True)

        n_cols = 5
        for i, g in enumerate(VARIABLE_GROUPS):
            v = BooleanVar(value=False)
            self.group_vars[g["group_key"]] = v
            row = i // n_cols
            col = i % n_cols
            cell = ttk.LabelFrame(grid, text="", padding=(3, 1))
            cell.grid(row=row, column=col, sticky="nsew", padx=2, pady=1)
            grid.rowconfigure(row, weight=1)
            grid.columnconfigure(col, weight=1)
            cb = ttk.Checkbutton(cell, text=f"{g['group_name']} ({len(g['variables'])}项)",
                                 variable=v)
            cb.pack(anchor="w")
            Label(cell, text=g["description"], font=("微软雅黑", 7), fg="#aaa",
                  wraplength=250, justify="left").pack(anchor="w", padx=(18, 0))

    def _build_custom_tab(self, parent):
        Label(parent, text="手动输入NHANES变量名（多个用逗号分隔）",
              font=("微软雅黑", 9)).pack(anchor="w", pady=5)
        Label(parent, text="格式: 直接输入变量名(如 RIAGENDR,BMXBMI,LBXTC)，或 表名:变量名(如 DEMO:RIDAGEYR)",
              font=("微软雅黑", 8), fg="#888").pack(anchor="w")

        self.custom_text = Text(parent, height=5, font=("Consolas", 10))
        self.custom_text.pack(fill="x", pady=5, padx=5)
        self.custom_text.insert("1.0", "RIAGENDR, RIDAGEYR, BMXBMI, LBXTC, LBXGLU")

        Label(parent, text="提示: 变量不存在时会自动跳过，不影响其他数据",
              font=("微软雅黑", 8), fg="#aaa").pack(anchor="w", padx=5)

    def _build_browser_tab(self, parent):
        """变量浏览器：选择周期后列出所有可用变量"""
        Label(parent, text="选择周期，列出该周期所有数据表的可用变量",
              font=("微软雅黑", 9)).pack(anchor="w", pady=2)

        row = ttk.Frame(parent)
        row.pack(fill="x", pady=3)
        Label(row, text="选择周期:").pack(side="left")
        self.browser_cycle_var = StringVar()
        cycle_options = [f"{c[0]} ({c[2]})" for c in NHANES_CYCLES]
        om = ttk.Combobox(row, textvariable=self.browser_cycle_var,
                          values=cycle_options, width=18, state="readonly")
        om.pack(side="left", padx=5)
        ttk.Button(row, text="获取变量", command=self._browser_fetch).pack(side="left")

        Label(parent, text="常见数据表:").pack(anchor="w", padx=5, pady=(5, 0))
        common_tables = "DEMO, BMX, BPX, BPQ, SMQ, ALQ, DIQ, PAQ, SLQ, HUQ, DPQ, RXQ_RX, " \
                        "TCHOL, GLU, INS, THYROD, BIOPRO, VID, HSCRP, CBC, PBCD, SEXHRM, DR1TOT"
        Label(parent, text=common_tables, font=("Consolas", 8), fg="#666",
              wraplength=700, justify="left").pack(anchor="w", padx=5)

        # 结果显示
        self.browser_result = Text(parent, height=5, font=("Consolas", 9))
        self.browser_result.pack(fill="both", expand=True, pady=5, padx=5)

    def _browser_fetch(self):
        """变量浏览器取数"""
        sel = self.browser_cycle_var.get()
        if not sel:
            messagebox.showwarning("提示", "请先选择一个周期")
            return
        suffix = sel.split("(")[1].split(")")[0]

        cycle_info = None
        for c in NHANES_CYCLES:
            if c[2] == suffix:
                cycle_info = c
                break
        if not cycle_info:
            return

        self.browser_result.delete("1.0", END)
        self.browser_result.insert(END, f"正在获取 {cycle_info[0]} 数据...\n")
        self.root.update()

        tables = ["DEMO", "BMX", "BPX", "BPQ", "SMQ", "ALQ", "DIQ", "PAQ", "SLQ",
                  "HUQ", "DPQ", "TCHOL", "GLU", "INS", "THYROD", "BIOPRO", "VID",
                  "HSCRP", "CBC", "PBCD", "SEXHRM", "DR1TOT", "RXQ_RX"]

        output = []
        for tn in tables:
            cols = self.engine.explorer.get_variables(tn, suffix, cycle_info[1])
            if cols:
                output.append(f"\n📄 {tn}_{suffix}.XPT ({len(cols)} 个变量):")
                # 每行5个
                for i in range(0, len(cols), 5):
                    chunk = cols[i:i+5]
                    output.append("  " + ", ".join(chunk))
            else:
                output.append(f"\n⏭️  {tn}_{suffix}: 无数据")

        self.browser_result.delete("1.0", END)
        self.browser_result.insert(END, f"=== {cycle_info[0]} 可用变量 ===\n")
        self.browser_result.insert(END, "\n".join(output))

    def _build_profile_tab(self, parent):
        """研究方案可视化编辑面板"""
        # 使用PanedWindow分上下区
        pw = ttk.PanedWindow(parent, orient="vertical")
        pw.pack(fill="both", expand=True)

        # === 上区：排除条件 ===
        top = ttk.LabelFrame(pw, text="排除条件设置", padding="5")
        pw.add(top, weight=1)

        # 疾病排除
        Label(top, text="排除疾病史（勾选=有该病史的排除）",
              font=("微软雅黑", 9, "bold")).pack(anchor="w")
        disease_frame = ttk.Frame(top)
        disease_frame.pack(fill="x", pady=2)
        self.profile_disease_vars = {}
        disease_list = sorted(DISEASE_EXCLUSION_CATEGORIES.keys())
        col1 = ttk.Frame(disease_frame)
        col2 = ttk.Frame(disease_frame)
        col1.pack(side="left", fill="x", expand=True)
        col2.pack(side="left", fill="x", expand=True)
        for i, dname in enumerate(disease_list):
            v = BooleanVar(value=False)
            self.profile_disease_vars[dname] = v
            target = col1 if i % 2 == 0 else col2
            ttk.Checkbutton(target, text=dname, variable=v).pack(anchor="w")

        ttk.Separator(top, orient="horizontal").pack(fill="x", pady=4)

        # 药物排除
        Label(top, text="排除用药（勾选=用该类药的排除）",
              font=("微软雅黑", 9, "bold")).pack(anchor="w")
        drug_frame = ttk.Frame(top)
        drug_frame.pack(fill="x", pady=2)
        self.profile_drug_vars = {}
        drug_list = sorted(DRUG_CLASSIFICATION.keys())
        col1 = ttk.Frame(drug_frame)
        col2 = ttk.Frame(drug_frame)
        col1.pack(side="left", fill="x", expand=True)
        col2.pack(side="left", fill="x", expand=True)
        for i, drname in enumerate(drug_list):
            v = BooleanVar(value=False)
            self.profile_drug_vars[drname] = v
            target = col1 if i % 2 == 0 else col2
            ttk.Checkbutton(target, text=drname, variable=v).pack(anchor="w")

        ttk.Separator(top, orient="horizontal").pack(fill="x", pady=4)

        # 年龄范围 + 范围排除快捷输入
        range_f = ttk.Frame(top)
        range_f.pack(fill="x")
        Label(range_f, text="年龄范围:", font=("微软雅黑", 9)).pack(side="left")
        self.profile_age_min = Entry(range_f, width=6)
        self.profile_age_min.insert(0, "18")
        self.profile_age_min.pack(side="left", padx=2)
        Label(range_f, text="~").pack(side="left")
        self.profile_age_max = Entry(range_f, width=6)
        self.profile_age_max.insert(0, "80")
        self.profile_age_max.pack(side="left", padx=2)
        Label(range_f, text="岁").pack(side="left")
        Label(range_f, text="  |  ", font=("微软雅黑", 8)).pack(side="left")

        # 衍生指标
        Label(range_f, text="衍生指标:", font=("微软雅黑", 9, "bold")).pack(side="left")
        idx_options = sorted(DERIVED_FORMULAS.keys())
        self.profile_idx_var = StringVar(value="无")
        idx_menu = ttk.Combobox(range_f, textvariable=self.profile_idx_var,
                                 values=["无"] + idx_options, width=16, state="readonly")
        idx_menu.pack(side="left", padx=2)
        self.profile_idx2_var = StringVar(value="无")
        idx_menu2 = ttk.Combobox(range_f, textvariable=self.profile_idx2_var,
                                  values=["无"] + idx_options, width=16, state="readonly")
        idx_menu2.pack(side="left", padx=2)

        ttk.Separator(top, orient="horizontal").pack(fill="x", pady=4)

        # 复合条件编辑器
        Label(top, text="复合筛选条件（支持 AND/OR 多组）",
              font=("微软雅黑", 9, "bold")).pack(anchor="w")
        Label(top, text="例: ALT(10~20) 且 AST(3~50) 为一组；组间可选 AND 或 OR 关系",
              font=("微软雅黑", 7), fg="#888").pack(anchor="w")

        cond_ctrl = ttk.Frame(top)
        cond_ctrl.pack(fill="x", pady=2)
        self.cond_between_logic = StringVar(value="AND")
        ttk.Radiobutton(cond_ctrl, text="组间 AND（满足全部组）", variable=self.cond_between_logic,
                        value="AND").pack(side="left", padx=2)
        ttk.Radiobutton(cond_ctrl, text="组间 OR（满足任意组）", variable=self.cond_between_logic,
                        value="OR").pack(side="left", padx=2)

        # 条件组列表容器
        cond_scroll = ttk.Frame(top)
        cond_scroll.pack(fill="both", expand=True)
        self.cond_canvas = Canvas(cond_scroll, height=120, highlightthickness=0)
        cond_sb = ttk.Scrollbar(cond_scroll, orient="vertical", command=self.cond_canvas.yview)
        self.cond_inner = ttk.Frame(self.cond_canvas)
        self.cond_inner.bind("<Configure>", lambda e: self.cond_canvas.configure(
            scrollregion=self.cond_canvas.bbox("all")))
        self.cond_canvas.create_window((0, 0), window=self.cond_inner, anchor="nw")
        self.cond_canvas.configure(yscrollcommand=cond_sb.set)
        self.cond_canvas.pack(side="left", fill="both", expand=True)
        cond_sb.pack(side="right", fill="y")
        _bind_mousewheel(self.cond_canvas)

        self.condition_rows = []  # 每行: [group_frame, logic_label, group_key, cond_widgets]
        ttk.Button(top, text="+ 添加条件组", command=self._add_condition_group).pack(anchor="w", pady=2)
        self._add_condition_group()  # 默认一组

        # === 下区：方案操作 ===
        bottom = ttk.LabelFrame(pw, text="方案操作", padding="5")
        pw.add(bottom, weight=0)

        btn_f = ttk.Frame(bottom)
        btn_f.pack(fill="x", pady=3)

        # 预置方案
        Label(btn_f, text="加载预置:", font=("微软雅黑", 9)).pack(side="left")
        self.profile_template_var = StringVar(value="无")
        tmpl_options = ["无"] + list(STUDY_TEMPLATES.keys())
        om = ttk.Combobox(btn_f, textvariable=self.profile_template_var,
                          values=tmpl_options, width=18, state="readonly")
        om.pack(side="left", padx=2)
        ttk.Button(btn_f, text="加载", width=5,
                   command=self._load_profile_template).pack(side="left", padx=2)

        ttk.Separator(btn_f, orient="vertical").pack(side="left", fill="y", padx=8)

        ttk.Button(btn_f, text="📥 保存方案", width=10,
                   command=self._save_visual_profile).pack(side="left", padx=2)
        ttk.Button(btn_f, text="📤 加载方案", width=10,
                   command=self._load_visual_profile).pack(side="left", padx=2)
        ttk.Button(btn_f, text="▶ 快速运行此方案", width=14,
                   command=self._run_visual_profile).pack(side="right", padx=2)

        self.profile_status = Label(bottom, text="就绪", font=("微软雅黑", 8), fg="#888")
        self.profile_status.pack(anchor="w", pady=2)

    def _load_profile_template(self):
        """从预置模板加载方案到可视化面板"""
        tname = self.profile_template_var.get()
        if tname == "无" or tname not in STUDY_TEMPLATES:
            return
        self._loaded_template_name = tname  # 记录模板名, 用于range_filters
        cfg = STUDY_TEMPLATES[tname]
        # 疾病排除
        for dname, v in self.profile_disease_vars.items():
            v.set(dname in cfg["exclusion"].get("disease_classes", []))
        # 药物排除
        for drname, v in self.profile_drug_vars.items():
            v.set(drname in cfg["exclusion"].get("medication_classes", []))
        # 年龄范围
        for criterion in cfg.get("inclusion", []):
            if "年龄" in criterion.get("var", ""):
                vmin = criterion.get("min", "")
                vmax = criterion.get("max", "")
                if vmin:
                    self.profile_age_min.delete(0, END)
                    self.profile_age_min.insert(0, str(vmin))
                if vmax:
                    self.profile_age_max.delete(0, END)
                    self.profile_age_max.insert(0, str(vmax))
        # 衍生指标
        idx_list = cfg.get("derived_indices", [])
        self.profile_idx_var.set(idx_list[0] if len(idx_list) > 0 else "无")
        self.profile_idx2_var.set(idx_list[1] if len(idx_list) > 1 else "无")
        self.profile_status.config(text=f"已加载模板: {tname}", fg="#2e7d32")
        # 复合条件
        self._load_compound_config_to_ui(cfg.get("compound_filters", {}))
        # 弹框确认
        detail_parts = []
        if cfg.get("inclusion"):
            detail_parts.append("纳入: " + "; ".join([f"{c['var']} {c.get('min','')}-{c.get('max','')}" for c in cfg["inclusion"]]))
        rf = cfg.get("exclusion", {}).get("range_filters", [])
        if rf:
            detail_parts.append(f"范围排除: {len(rf)}项")
        if cfg.get("exclusion", {}).get("disease_classes"):
            detail_parts.append("疾病排除: " + ", ".join(cfg["exclusion"]["disease_classes"]))
        if cfg.get("exclusion", {}).get("medication_classes"):
            detail_parts.append("药物排除: " + ", ".join(cfg["exclusion"]["medication_classes"]))
        detail = "\n".join(detail_parts) if detail_parts else "无特殊筛选条件"
        messagebox.showinfo("模板已加载", f"已加载模板: {tname}\n\n{detail}")

    def _build_profile_from_ui(self):
        """从GUI可视化面板构建方案配置"""
        inclusion = []
        try:
            age_min = int(self.profile_age_min.get())
            age_max = int(self.profile_age_max.get())
            inclusion.append({"var": "年龄(Age,岁)", "min": age_min, "max": age_max})
        except:
            pass

        disease_classes = [d for d, v in self.profile_disease_vars.items() if v.get()]
        medication_classes = [d for d, v in self.profile_drug_vars.items() if v.get()]

        # 如果UI复选框为空但模板有值，从模板读取
        if (not disease_classes or not medication_classes) and self._loaded_template_name and self._loaded_template_name in STUDY_TEMPLATES:
            tmpl_exc = STUDY_TEMPLATES[self._loaded_template_name].get("exclusion", {})
            if not disease_classes:
                disease_classes = tmpl_exc.get("disease_classes", [])
            if not medication_classes:
                medication_classes = tmpl_exc.get("medication_classes", [])

        derived = []
        for idx_var in [self.profile_idx_var, self.profile_idx2_var]:
            val = idx_var.get()
            if val and val != "无":
                derived.append(val)

        # 从已加载模板读取range_filters（UI没有对应控件，直接从STUDY_TEMPLATES拿）
        range_filters = []
        if self._loaded_template_name and self._loaded_template_name in STUDY_TEMPLATES:
            range_filters = STUDY_TEMPLATES[self._loaded_template_name].get("exclusion", {}).get("range_filters", [])
            # 同时补全inclusion中的BMI等条件
            tmpl_inclusion = STUDY_TEMPLATES[self._loaded_template_name].get("inclusion", [])
            existing_vars = {c["var"] for c in inclusion}
            for tc in tmpl_inclusion:
                if tc["var"] not in existing_vars:
                    inclusion.append(tc)

        return {
            "inclusion": inclusion,
            "exclusion": {
                "medication_classes": medication_classes,
                "disease_classes": disease_classes,
                "range_filters": range_filters,  # 现在有了
                "drug_audit": True,
            },
            "compound_filters": self._build_compound_config_from_ui(),
            "derived_indices": derived,
            "disease_classes": disease_classes,
            "medication_classes": medication_classes,
        }

    def _save_visual_profile(self):
        """保存可视化方案为JSON"""
        cfg = self._build_profile_from_ui()
        path = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON", "*.json")],
            initialdir=PROFILES_DIR,
            initialfile="my_study_profile.json"
        )
        if not path:
            return
        # 构建完整方案
        profile = StudyProfile()
        selected_cycles = [s for s, v in self.cycle_vars.items() if v.get()]
        selected_groups = [k for k, v in self.group_vars.items() if v.get()]
        profile.config["study_name"] = os.path.splitext(os.path.basename(path))[0]
        profile.config["data"]["cycles"] = selected_cycles
        profile.config["data"]["groups"] = selected_groups
        profile.config["data"]["convert_units"] = self.convert_var.get()
        profile.config["data"]["auto_aggregate"] = self.aggregate_var.get()
        profile.config["data"]["include_mortality"] = self.mort_var.get()
        profile.config["inclusion"] = cfg["inclusion"]
        profile.config["exclusion"]["disease_classes"] = cfg["disease_classes"]
        profile.config["exclusion"]["medication_classes"] = cfg["medication_classes"]
        profile.config["compound_filters"] = cfg.get("compound_filters", {})
        profile.config["derived_indices"] = cfg["derived_indices"]
        if profile.save(path):
            self.profile_status.config(text=f"已保存: {os.path.basename(path)}", fg="#1565c0")
            self._log(f"📋 研究方案已保存: {path}")
        else:
            self.profile_status.config(text="保存失败", fg="#c62828")

    # ────────────────────────────────────────────────────────────
    # 模板管理方法（右侧控制面板的"方案模板"区域）
    # ────────────────────────────────────────────────────────────

    def _refresh_custom_templates(self):
        """刷新自定义模板列表"""
        self.custom_templates = []
        profiles = StudyProfile.list_profiles()
        for p in profiles:
            self.custom_templates.append({"name": p["name"], "path": p["path"], "desc": p["desc"]})

    def _load_template(self):
        """加载选中的模板到当前UI设置"""
        sel = self.tmpl_var.get()
        if not sel or sel == "— 选择模板 —":
            return

        # 检查是预置模板还是自定义模板
        if sel in STUDY_TEMPLATES:
            tmpl_cfg = STUDY_TEMPLATES[sel]
            self._log(f"📋 加载预置模板: {sel}")
            self._loaded_template_name = sel
        else:
            # 自定义模板，从JSON加载
            for t in self.custom_templates:
                if t["name"] == sel:
                    profile = StudyProfile(t["path"])
                    if profile.load():
                        tmpl_cfg = profile.config
                        self._log(f"📋 加载自定义模板: {sel}")
                        self._loaded_template_name = None
                        break
            else:
                self._log(f"⚠️ 未找到模板: {sel}")
                return

        # 应用到UI：周期
        cycles = tmpl_cfg.get("data", {}).get("cycles", [])
        if cycles:
            for c_s, v in self.cycle_vars.items():
                v.set(c_s in cycles)

        # 应用到UI：变量组
        groups = tmpl_cfg.get("data", {}).get("groups", [])
        if groups:
            for g_k, v in self.group_vars.items():
                v.set(g_k in groups)

        # 应用到UI：选项
        if "convert_units" in tmpl_cfg.get("data", {}):
            self.convert_var.set(tmpl_cfg["data"]["convert_units"])
        if "auto_aggregate" in tmpl_cfg.get("data", {}):
            self.aggregate_var.set(tmpl_cfg["data"]["auto_aggregate"])
        if tmpl_cfg.get("data", {}).get("include_mortality", False):
            self.mort_var.set(True)

        # 应用到UI：疾病排除、药物排除、年龄范围
        exc = tmpl_cfg.get("exclusion", {})
        if hasattr(self, 'profile_disease_vars'):
            for dname, v in self.profile_disease_vars.items():
                v.set(dname in exc.get("disease_classes", []))
        if hasattr(self, 'profile_drug_vars'):
            for drname, v in self.profile_drug_vars.items():
                v.set(drname in exc.get("medication_classes", []))
        if hasattr(self, 'profile_age_min') and hasattr(self, 'profile_age_max'):
            for criterion in tmpl_cfg.get("inclusion", []):
                if "年龄" in criterion.get("var", ""):
                    vmin = criterion.get("min", "")
                    vmax = criterion.get("max", "")
                    if vmin:
                        self.profile_age_min.delete(0, END)
                        self.profile_age_min.insert(0, str(vmin))
                    if vmax:
                        self.profile_age_max.delete(0, END)
                        self.profile_age_max.insert(0, str(vmax))
        # 应用到UI：衍生指标
        idx_list = tmpl_cfg.get("derived_indices", [])
        if hasattr(self, 'profile_idx_var'):
            self.profile_idx_var.set(idx_list[0] if len(idx_list) > 0 else "无")
            self.profile_idx2_var.set(idx_list[1] if len(idx_list) > 1 else "无")

        self._log(f"   ✅ 已应用: 周期={cycles}, 变量组={groups}")

        # 弹框确认
        detail_parts = []
        inc = tmpl_cfg.get("inclusion", [])
        if inc:
            detail_parts.append("纳入: " + "; ".join([f"{c['var']} {c.get('min','')}-{c.get('max','')}" for c in inc]))
        exc = tmpl_cfg.get("exclusion", {})
        rf = exc.get("range_filters", [])
        if rf:
            detail_parts.append(f"范围排除: {len(rf)}项")
        if exc.get("disease_classes"):
            detail_parts.append("疾病排除: " + ", ".join(exc["disease_classes"]))
        if exc.get("medication_classes"):
            detail_parts.append("药物排除: " + ", ".join(exc["medication_classes"]))
        detail = "\n".join(detail_parts) if detail_parts else "无特殊筛选条件"
        messagebox.showinfo("模板已加载", f"已加载: {sel}\n\n{detail}")

    def _save_as_template(self):
        """将当前UI设置保存为自定义模板"""
        selected_cycles = [s for s, v in self.cycle_vars.items() if v.get()]
        selected_groups = [k for k, v in self.group_vars.items() if v.get()]

        name = simpledialog.askstring("保存模板", "请输入模板名称:")
        if not name:
            return
        name = name.strip()
        if not name:
            return

        cfg = {
            "study_name": name,
            "description": f"自定义模板 - {name}",
            "data": {
                "cycles": selected_cycles,
                "groups": selected_groups,
                "convert_units": self.convert_var.get(),
                "auto_aggregate": self.aggregate_var.get(),
                "include_mortality": self.mort_var.get(),
            },
            "inclusion": [],
            "exclusion": {"medication_classes": [], "range_filters": [], "disease_classes": [], "drug_audit": True},
            "derived_indices": [],
        }

        path = os.path.join(PROFILES_DIR, f"{name}.json")
        with open(path, "w", encoding="utf-8") as f:
            json.dump(cfg, f, ensure_ascii=False, indent=2)

        self._refresh_custom_templates()
        tmpl_options = ["— 选择模板 —"] + list(STUDY_TEMPLATES.keys()) + [t["name"] for t in self.custom_templates]
        for w in self.root.winfo_children():
            for child in w.winfo_children():
                for c2 in child.winfo_children():
                    if isinstance(c2, ttk.Combobox) and c2.cget("values") and "选择模板" in str(c2.cget("values")[0]):
                        c2.configure(values=tmpl_options)
                        break
        self.tmpl_var.set(name)
        self._log(f"📋 模板已保存: {name}")

    def _manage_templates(self):
        """管理自定义模板：编辑名称、删除"""
        if not self.custom_templates:
            messagebox.showinfo("提示", "暂无自定义模板")
            return

        # 弹出选择对话框
        dialog = tk.Toplevel(self.root)
        dialog.title("管理自定义模板")
        dialog.geometry("450x350")
        dialog.transient(self.root)
        dialog.grab_set()

        Label(dialog, text="自定义模板列表", font=("微软雅黑", 11, "bold")).pack(pady=8)

        frame = ttk.Frame(dialog)
        frame.pack(fill="both", expand=True, padx=10, pady=5)

        # 列表框
        lb = tk.Listbox(frame, font=("微软雅黑", 10))
        sb = ttk.Scrollbar(frame, orient="vertical", command=lb.yview)
        lb.configure(yscrollcommand=sb.set)
        lb.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")

        for t in self.custom_templates:
            lb.insert(tk.END, t["name"])

        def delete_selected():
            sel = lb.curselection()
            if not sel:
                return
            idx = sel[0]
            name = self.custom_templates[idx]["name"]
            path = self.custom_templates[idx]["path"]
            if messagebox.askyesno("确认删除", f"确定删除模板「{name}」？"):
                try:
                    os.remove(path)
                    lb.delete(idx)
                    self.custom_templates.pop(idx)
                    self._log(f"🗑 模板已删除: {name}")
                    # 更新下拉
                    tmpl_options = ["— 选择模板 —"] + list(STUDY_TEMPLATES.keys()) + [t["name"] for t in self.custom_templates]
                    for w in self.root.winfo_children():
                        for child in w.winfo_children():
                            for c2 in child.winfo_children():
                                if isinstance(c2, ttk.Combobox) and "选择模板" in str(c2.cget("values")[0]):
                                    c2.configure(values=tmpl_options)
                                    break
                except Exception as e:
                    messagebox.showerror("错误", f"删除失败: {e}")

        def rename_selected():
            sel = lb.curselection()
            if not sel:
                return
            idx = sel[0]
            old_name = self.custom_templates[idx]["name"]
            old_path = self.custom_templates[idx]["path"]
            new_name = simpledialog.askstring("重命名", f"将「{old_name}」重命名为:")
            if not new_name:
                return
            new_name = new_name.strip()
            if not new_name:
                return
            try:
                # 更新JSON文件
                with open(old_path, "r", encoding="utf-8") as f:
                    cfg = json.load(f)
                cfg["study_name"] = new_name
                new_path = os.path.join(PROFILES_DIR, f"{new_name}.json")
                with open(new_path, "w", encoding="utf-8") as f:
                    json.dump(cfg, f, ensure_ascii=False, indent=2)
                os.remove(old_path)
                self.custom_templates[idx]["name"] = new_name
                self.custom_templates[idx]["path"] = new_path
                lb.delete(idx)
                lb.insert(idx, new_name)
                self._log(f"✏️ 模板已重命名: {old_name} -> {new_name}")
                # 更新下拉
                tmpl_options = ["— 选择模板 —"] + list(STUDY_TEMPLATES.keys()) + [t["name"] for t in self.custom_templates]
                for w in self.root.winfo_children():
                    for child in w.winfo_children():
                        for c2 in child.winfo_children():
                            if isinstance(c2, ttk.Combobox) and "选择模板" in str(c2.cget("values")[0]):
                                c2.configure(values=tmpl_options)
                                break
            except Exception as e:
                messagebox.showerror("错误", f"重命名失败: {e}")

        btn_f = ttk.Frame(dialog)
        btn_f.pack(fill="x", pady=8)
        ttk.Button(btn_f, text="✏️ 重命名", command=rename_selected).pack(side="left", padx=5)
        ttk.Button(btn_f, text="🗑 删除", command=delete_selected).pack(side="left", padx=5)
        ttk.Button(btn_f, text="关闭", command=dialog.destroy).pack(side="right", padx=5)

    def _load_visual_profile(self):
        """加载JSON方案到可视化面板"""
        path = filedialog.askopenfilename(
            filetypes=[("JSON", "*.json")],
            initialdir=PROFILES_DIR
        )
        if not path:
            return
        profile = StudyProfile()
        if not profile.load(path):
            self.profile_status.config(text="加载失败", fg="#c62828")
            return
        cfg = profile.config
        # 疾病排除
        for dname, v in self.profile_disease_vars.items():
            v.set(dname in cfg["exclusion"].get("disease_classes", []))
        # 药物排除
        for drname, v in self.profile_drug_vars.items():
            v.set(drname in cfg["exclusion"].get("medication_classes", []))
        # 年龄
        for criterion in cfg.get("inclusion", []):
            if "年龄" in criterion.get("var", ""):
                vmin = criterion.get("min", "")
                vmax = criterion.get("max", "")
                if vmin:
                    self.profile_age_min.delete(0, END)
                    self.profile_age_min.insert(0, str(vmin))
                if vmax:
                    self.profile_age_max.delete(0, END)
                    self.profile_age_max.insert(0, str(vmax))
        # 衍生指标
        idx_list = cfg.get("derived_indices", [])
        self.profile_idx_var.set(idx_list[0] if len(idx_list) > 0 else "无")
        self.profile_idx2_var.set(idx_list[1] if len(idx_list) > 1 else "无")
        # 复合条件
        self._load_compound_config_to_ui(cfg.get("compound_filters", {}))
        # 同步勾选周期和变量组
        for s, v in self.cycle_vars.items():
            v.set(s in cfg["data"].get("cycles", []))
        for k, v in self.group_vars.items():
            v.set(k in cfg["data"].get("groups", []))
        self.convert_var.set(cfg["data"].get("convert_units", True))
        self.aggregate_var.set(cfg["data"].get("auto_aggregate", True))
        self.mort_var.set(cfg["data"].get("include_mortality", False))
        self.profile_status.config(text=f"已加载: {os.path.basename(path)}", fg="#2e7d32")
        self._log(f"📋 研究方案已加载: {path}")

    def _run_visual_profile(self):
        """快速运行当前可视化方案"""
        cfg = self._build_profile_from_ui()
        selected_cycles = [s for s, v in self.cycle_vars.items() if v.get()]
        selected_groups = [k for k, v in self.group_vars.items() if v.get()]
        if not selected_cycles:
            messagebox.showwarning("提示", "请先选择周期（左侧）")
            return
        if not selected_groups:
            messagebox.showwarning("提示", "请先选择变量组（预定义变量组Tab）")
            return

        profile = StudyProfile()
        profile.config["study_name"] = "可视化方案"
        profile.config["data"]["cycles"] = selected_cycles
        profile.config["data"]["groups"] = selected_groups
        profile.config["data"]["convert_units"] = self.convert_var.get()
        profile.config["data"]["auto_aggregate"] = self.aggregate_var.get()
        profile.config["data"]["include_mortality"] = self.mort_var.get()
        profile.config["inclusion"] = cfg["inclusion"]
        profile.config["exclusion"]["disease_classes"] = cfg["disease_classes"]
        profile.config["exclusion"]["medication_classes"] = cfg["medication_classes"]
        profile.config["derived_indices"] = cfg["derived_indices"]

        out_path = self.path_var.get().strip()
        if not out_path:
            out_path = filedialog.askdirectory(
                initialdir=self.path_var.get() or os.path.expanduser("~"),
                title="选择保存文件夹（自动按年份建子目录）")
            if not out_path:
                return
            self.path_var.set(out_path)

        msg = f"周期: {len(selected_cycles)}个\n"
        msg += f"排除疾病: {cfg['disease_classes'] or '无'}\n"
        msg += f"排除药物: {cfg['medication_classes'] or '无'}\n"
        cf = cfg.get("compound_filters", {})
        if cf and cf.get("groups"):
            n_total_conds = sum(len(g.get("conditions",[])) for g in cf["groups"])
            msg += f"复合条件: {len(cf['groups'])}组共{n_total_conds}项\n"
        msg += f"衍生指标: {cfg['derived_indices'] or '无'}\n"
        msg += f"输出: {out_path}\n\n确定运行？"
        if not messagebox.askyesno("运行方案", msg):
            return

        self.start_btn.configure(state="disabled")
        self._log("="*50)
        self._log(f"📋 运行研究方案: 可视化方案")

        t = threading.Thread(
            target=self._run_task,
            args=(selected_cycles, selected_groups, "", self.mort_var.get(),
                  out_path, self.convert_var.get(), self.aggregate_var.get(),
                  profile),
            daemon=True
        )
        t.start()
        self._check_status()

    # ========== 复合条件编辑器辅助方法 ==========

    def _add_condition_group(self, within_logic="AND", conditions=None):
        """添加一个条件组到UI"""
        gidx = len(self.condition_rows)
        gf = ttk.LabelFrame(self.cond_inner, text=f"组 {gidx+1}", padding="3")
        gf.pack(fill="x", pady=2, padx=2)

        row0 = ttk.Frame(gf)
        row0.pack(fill="x")
        logic_var = StringVar(value=within_logic)
        ttk.Radiobutton(row0, text="组内 AND（全部满足）", variable=logic_var,
                        value="AND").pack(side="left", padx=1)
        ttk.Radiobutton(row0, text="组内 OR（任一满足）", variable=logic_var,
                        value="OR").pack(side="left", padx=1)
        # 使用闭包绑定当前索引
        def make_del_cb(idx):
            return lambda: [self.condition_rows.pop(idx), gf.destroy()] if idx < len(self.condition_rows) else None
        ttk.Button(row0, text="✕ 删除本组", width=10,
                   command=make_del_cb(gidx)).pack(side="right", padx=2)

        cond_list = []  # [(cond_frame, var_menu, op_menu, v1_entry, v2_entry)]
        cond_container = ttk.Frame(gf)
        cond_container.pack(fill="x")

        def add_cond(cfg=None):
            cf = ttk.Frame(cond_container)
            cf.pack(fill="x", pady=1)
            # 变量名输入
            var_entry = ttk.Combobox(cf, width=20, values=sorted(list(set(
                v["var_label"] for g in VARIABLE_GROUPS for v in g["variables"] if not v.get("is_id")))))
            var_entry.set(cfg.get("var", "BMI(kg/m²)") if cfg else "")
            var_entry.pack(side="left", padx=1)
            # 操作符
            op_var = StringVar(value=cfg.get("op", "区间") if cfg else "区间")
            op_menu = ttk.Combobox(cf, textvariable=op_var, values=FILTER_OPERATORS, width=5, state="readonly")
            op_menu.pack(side="left", padx=1)
            # 值1
            v1 = Entry(cf, width=7)
            v1.insert(0, str(cfg.get("min", cfg.get("value", ""))) if cfg else "")
            v1.pack(side="left", padx=1)
            # 值2（仅区间用）
            Label(cf, text="~").pack(side="left")
            v2 = Entry(cf, width=7)
            v2.insert(0, str(cfg.get("max", "")) if cfg else "")
            v2.pack(side="left", padx=1)
            # 删除此条件
            ttk.Button(cf, text="✕", width=2, command=cf.destroy).pack(side="right", padx=2)
            cond_list.append((cf, var_entry, op_var, v1, v2))

        # 默认加两条空条件
        if conditions:
            for c in conditions:
                add_cond(c)
        else:
            add_cond({"var": "", "op": "区间"})
            add_cond({"var": "", "op": "区间"})

        ttk.Button(gf, text="+ 添加条件", command=add_cond).pack(anchor="w", pady=1)

        self.condition_rows.append([gf, logic_var, cond_list])

    def _build_compound_config_from_ui(self):
        """从UI构建复合条件配置"""
        groups = []
        for gf, logic_var, cond_list in self.condition_rows:
            if not gf.winfo_exists():
                continue
            conditions = []
            for cf, var_entry, op_var, v1, v2 in cond_list:
                if not cf.winfo_exists():
                    continue
                var_name = var_entry.get().strip()
                op = op_var.get()
                if not var_name:
                    continue
                cond = {"var": var_name, "op": op}
                if op == "区间":
                    try:
                        cond["min"] = float(v1.get())
                        cond["max"] = float(v2.get())
                    except:
                        continue
                else:
                    try:
                        cond["value"] = float(v1.get())
                    except:
                        continue
                conditions.append(cond)
            if conditions:
                groups.append({"within_group": logic_var.get(), "conditions": conditions})
        return {"between_groups": self.cond_between_logic.get(), "groups": groups} if groups else {}

    def _load_compound_config_to_ui(self, filter_config):
        """加载复合条件配置到UI"""
        if not filter_config or "groups" not in filter_config:
            return
        # 清除已有组
        for gf, _, _ in self.condition_rows:
            try:
                gf.destroy()
            except:
                pass
        self.condition_rows = []
        self.cond_between_logic.set(filter_config.get("between_groups", "AND"))
        for g in filter_config["groups"]:
            self._add_condition_group(
                within_logic=g.get("within_group", "AND"),
                conditions=g.get("conditions", [])
            )
        if not self.condition_rows:
            self._add_condition_group()

    def _build_control_panel(self, parent):
        # 右侧面板内容太多，包装到可滚动画布中
        canvas = Canvas(parent, highlightthickness=0)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scroll_frame = ttk.Frame(canvas)
        scroll_frame.pack(fill="both", expand=True)
        scroll_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scroll_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        # 鼠标滚轮绑定
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        canvas.bind_all("<MouseWheel>", _on_mousewheel)
        parent = scroll_frame  # 所有子控件都放到scroll_frame内
        # 输出设置
        Label(parent, text="输出设置", font=("微软雅黑", 10, "bold")).pack(anchor="w", pady=2)

        ttk.Label(parent, text="格式:").pack(anchor="w")
        self.fmt_var = StringVar(value="CSV (.csv)")
        ttk.Combobox(parent, textvariable=self.fmt_var,
                     values=["CSV (.csv)", "Excel (.xlsx)"],
                     state="readonly", width=16).pack(anchor="w", pady=2)

        ttk.Label(parent, text="保存到文件夹:").pack(anchor="w")
        p_f = ttk.Frame(parent)
        p_f.pack(fill="x", pady=2)
        self.path_var = StringVar(
            value=os.path.join(os.path.expanduser("~"), "Desktop"))
        ttk.Entry(p_f, textvariable=self.path_var, width=18).pack(side="left", fill="x", expand=True)
        ttk.Button(p_f, text="浏览", width=5, command=self._browse_path).pack(side="right")

        ttk.Separator(parent, orient="horizontal").pack(fill="x", pady=8)

        # 额外选项
        Label(parent, text="额外选项", font=("微软雅黑", 10, "bold")).pack(anchor="w", pady=2)
        self.mort_var = BooleanVar(value=False)
        ttk.Checkbutton(parent, text="💀 包含死亡率数据(自动下载CDC链接文件)",
                        variable=self.mort_var).pack(anchor="w", pady=2)
        self.convert_var = BooleanVar(value=True)
        ttk.Checkbutton(parent, text="📐 单位换算(美制→中国常用单位)",
                        variable=self.convert_var).pack(anchor="w", pady=2)
        self.aggregate_var = BooleanVar(value=True)
        ttk.Checkbutton(parent, text="🔄 自动聚合(吃鱼/吸烟/饮酒等合并为综合指标)",
                        variable=self.aggregate_var).pack(anchor="w", pady=2)

        ttk.Separator(parent, orient="horizontal").pack(fill="x", pady=8)

        # 方案模板管理
        Label(parent, text="方案模板", font=("微软雅黑", 10, "bold")).pack(anchor="w", pady=2)
        tmpl_frame = ttk.Frame(parent)
        tmpl_frame.pack(fill="x")
        self.tmpl_var = StringVar(value="— 选择模板 —")
        tmpl_options = ["— 选择模板 —"] + list(STUDY_TEMPLATES.keys())
        # 加载自定义模板
        self._refresh_custom_templates()
        tmpl_options += [t["name"] for t in self.custom_templates]
        tmpl_om = ttk.Combobox(tmpl_frame, textvariable=self.tmpl_var,
                               values=tmpl_options, width=18, state="readonly")
        tmpl_om.pack(side="left", padx=1)
        ttk.Button(tmpl_frame, text="加载", width=5,
                   command=self._load_template).pack(side="left", padx=1)
        ttk.Button(tmpl_frame, text="💾 保存", width=5,
                   command=self._save_as_template).pack(side="left", padx=1)
        ttk.Button(tmpl_frame, text="⚙ 管理", width=5,
                   command=self._manage_templates).pack(side="left", padx=1)

        ttk.Separator(parent, orient="horizontal").pack(fill="x", pady=8)

        # 操作按钮
        self.start_btn = ttk.Button(parent, text="▶ 开始处理", command=self._start,
                                     width=18)
        self.start_btn.pack(fill="x", pady=3)

        self.preview_btn = ttk.Button(parent, text="👁 预览结果", command=self._preview,
                                       width=18, state="disabled")
        self.preview_btn.pack(fill="x", pady=3)

        ttk.Separator(parent, orient="horizontal").pack(fill="x", pady=5)
        Label(parent, text="🔍 数据验证", font=("微软雅黑", 9, "bold")).pack(anchor="w", pady=1)

        self.verify_strict_btn = ttk.Button(parent, text="✅ 0.1N严格验证",
                                             command=self._run_verify_strict,
                                             width=18, state="disabled")
        self.verify_strict_btn.pack(fill="x", pady=2)
        self.verify_lenient_btn = ttk.Button(parent, text="🔶 0.1N宽松验证",
                                              command=self._run_verify_lenient,
                                              width=18, state="disabled")
        self.verify_lenient_btn.pack(fill="x", pady=2)

        ttk.Separator(parent, orient="horizontal").pack(fill="x", pady=5)

        Label(parent, text="本地数据库", font=("微软雅黑", 10, "bold")).pack(anchor="w", pady=2)
        self.db_btn = ttk.Button(parent, text="🗄️ 构建本地数据库\n(一次性下载全部周期)",
                                  command=self._build_db, width=18)
        self.db_btn.pack(fill="x", pady=3)
        self.db_import_btn = ttk.Button(parent, text="📂 导入本地数据库\n(选择已有.db文件)",
                                        command=self._import_db, width=18)
        self.db_import_btn.pack(fill="x", pady=1)
        self.db_query_btn = ttk.Button(parent, text="🔍 从本地数据库查询",
                                       command=self._query_db, width=18,
                                       state="disabled")
        self.db_query_btn.pack(fill="x", pady=3)

        ttk.Separator(parent, orient="horizontal").pack(fill="x", pady=8)

        # 状态
        self.status_var = StringVar(value="就绪")
        ttk.Label(parent, textvariable=self.status_var,
                  font=("微软雅黑", 9)).pack(anchor="w")

    def _build_log_panel(self, parent):
        f = ttk.Frame(parent)
        f.pack(fill="both", expand=True)
        self.log_text = Text(f, height=7, font=("Consolas", 9),
                             bg="#1e1e1e", fg="#d4d4d4", wrap="word", state="disabled")
        self.log_text.pack(side="left", fill="both", expand=True)
        sb = ttk.Scrollbar(f, orient="vertical", command=self.log_text.yview)
        sb.pack(side="right", fill="y")
        self.log_text.configure(yscrollcommand=sb.set)

    def _log(self, msg):
        self.log_text.configure(state="normal")
        ts = datetime.now().strftime("%H:%M:%S")
        self.log_text.insert(END, f"[{ts}] {msg}\n")
        self.log_text.see(END)
        self.log_text.configure(state="disabled")
        self.root.update_idletasks()

    def _browse_path(self):
        p = filedialog.askdirectory(
            initialdir=self.path_var.get() or os.path.expanduser("~"),
            title="选择保存文件夹（自动按年份建子目录）")
        if p:
            self.path_var.set(p)

    def _select_sci(self):
        """推荐SCI常用变量"""
        sci = {"demo_core", "bmx", "bpx", "bpq", "smq", "alq", "diq", "paq",
               "slq", "dpq", "lipid", "glu", "thy", "kidney", "liver",
               "vit", "inf", "metal", "huq", "diet", "rxq"}
        for k, v in self.group_vars.items():
            v.set(k in sci)

    def _start(self):
        selected_cycles = [s for s, v in self.cycle_vars.items() if v.get()]
        selected_groups = [k for k, v in self.group_vars.items() if v.get()]
        custom_vars = self.custom_text.get("1.0", END).strip() if hasattr(self, 'custom_text') else ""
        include_mort = self.mort_var.get()
        convert_units = self.convert_var.get()
        auto_aggregate = self.aggregate_var.get()
        out_path = self.path_var.get().strip()

        if not selected_cycles:
            messagebox.showwarning("提示", "请选择至少一个调查周期！")
            return
        if not selected_groups and not custom_vars:
            messagebox.showwarning("提示", "请选择至少一个变量组，或输入自定义变量！")
            return
        if not out_path:
            messagebox.showwarning("提示", "请设置保存文件夹！")
            return
        # 确保是目录（去掉多余字符串）
        if not out_path.endswith((":", "\\", "/")) and not os.path.isdir(out_path):
            # 可能用户拖了个文件路径进来, 取其目录
            d = os.path.dirname(out_path)
            if d and os.path.isdir(d):
                out_path = d
            else:
                messagebox.showwarning("提示", "保存路径必须是文件夹！")
                return
        self.path_var.set(out_path)

        # 从可视化面板构建方案配置
        cfg = self._build_profile_from_ui()
        profile = StudyProfile()
        profile.config["study_name"] = "可视化方案"
        profile.config["data"]["cycles"] = selected_cycles
        profile.config["data"]["groups"] = selected_groups
        profile.config["data"]["convert_units"] = convert_units
        profile.config["data"]["auto_aggregate"] = auto_aggregate
        profile.config["data"]["include_mortality"] = include_mort
        profile.config["inclusion"] = cfg["inclusion"]
        profile.config["exclusion"]["disease_classes"] = cfg["disease_classes"]
        profile.config["exclusion"]["medication_classes"] = cfg["medication_classes"]
        profile.config["exclusion"]["range_filters"] = cfg["exclusion"]["range_filters"]
        profile.config["derived_indices"] = cfg["derived_indices"]
        profile.config["compound_filters"] = cfg.get("compound_filters", {})

        # 确认
        cycles_str = ", ".join([f"{c[0]}" for c in NHANES_CYCLES if c[2] in selected_cycles])
        groups_str = ", ".join([g["group_name"] for g in VARIABLE_GROUPS
                                 if g["group_key"] in selected_groups])
        msg = f"周期: {cycles_str}\n指标组: {groups_str or '(仅自定义变量)'}\n"
        if custom_vars.strip():
            msg += f"自定义: {custom_vars.strip()[:60]}...\n"
        if include_mort:
            msg += "含死亡率数据\n"
        if convert_units:
            msg += "含单位换算(美→中)\n"
        if auto_aggregate:
            msg += "含自动聚合(吃鱼/吸烟/饮酒等)\n"
        # 显示筛选条件
        if cfg["inclusion"]:
            msg += "纳入条件: " + "; ".join([f"{c['var']} {c.get('min','')}-{c.get('max','')}" for c in cfg["inclusion"]]) + "\n"
        if cfg["exclusion"]["range_filters"]:
            msg += f"范围排除: {len(cfg['exclusion']['range_filters'])}项\n"
        if cfg["disease_classes"]:
            msg += f"疾病排除: {', '.join(cfg['disease_classes'])}\n"
        if cfg["medication_classes"]:
            msg += f"药物排除: {', '.join(cfg['medication_classes'])}\n"
        msg += f"输出: {out_path}\n\n确定开始？"
        if not messagebox.askyesno("确认", msg):
            return

        self.start_btn.configure(state="disabled")
        self.progress.configure(value=0)
        self._log("="*50)
        self._log("🚀 任务开始...")

        t = threading.Thread(target=self._run_task,
                             args=(selected_cycles, selected_groups,
                                   custom_vars, include_mort, out_path,
                                   convert_units, auto_aggregate),
                             kwargs={"profile": profile},
                             daemon=True)
        t.start()
        self._check_status()

    def _run_task(self, cycles, groups, custom_vars, mort, out_path,
                  convert_units=False, auto_aggregate=False, profile=None):
        try:
            result = self.engine.run(cycles, groups, custom_vars, mort, out_path,
                                     convert_units, auto_aggregate, profile=profile)
            self.last_result = result
        except Exception as e:
            self.last_result = {"success": False, "error": str(e)}

    def _check_status(self):
        if self.engine.is_running:
            self.status_var.set("⏳ 处理中...")
            self.root.after(500, self._check_status)
        else:
            self.start_btn.configure(state="normal")
            r = self.last_result
            if r and r.get("success"):
                self.status_var.set(f"✅ 完成! {r['rows']}条 x {r['cols']}列")
                self.progress.configure(value=100)
                self.preview_btn.configure(state="normal")
                if _VALIDATION_AVAILABLE:
                    self.verify_strict_btn.configure(state="normal")
                    self.verify_lenient_btn.configure(state="normal")
                self._log(f"\n✅ 成功导出: {r['file_path']}")
                self._log(f"   共 {r['rows']} 条记录, {r['cols']} 列")
                # 弹摘要窗口
                self._show_completion_summary(r)
            else:
                self.status_var.set(f"❌ 失败")
                self._log(f"\n❌ 处理失败: {r.get('error', '未知错误') if r else '无返回'}")

    def _show_completion_summary(self, result):
        """运行完成后弹出摘要窗口"""
        if not result:
            return
        qc = result.get("qc_report", {}) or {}
        rows = result.get("rows", 0)
        cols = result.get("cols", 0)
        fpath = result.get("file_path", "")

        # 构建摘要文本
        lines = []
        lines.append(f"✅ NHANES 数据处理完成\n")
        lines.append(f"📊 总记录数: {rows:,} 人")
        lines.append(f"📋 总变量数: {cols} 列")
        lines.append(f"💾 输出文件: {os.path.basename(fpath) if fpath else '未导出'}")
        # 周期分布
        cycle_dist = qc.get("周期分布", {})
        if cycle_dist:
            lines.append(f"\n📅 周期分布:")
            for k, v in list(cycle_dist.items())[:6]:
                lines.append(f"   {k}: {v} 人")
            if len(cycle_dist) > 6:
                lines.append(f"   ... 共 {len(cycle_dist)} 个周期")
        # 高缺失率
        high_miss = {k: v for k, v in qc.get("变量缺失率", {}).items() if v > 50}
        if high_miss:
            lines.append(f"\n⚠️ 高缺失率变量(>50%): {len(high_miss)} 个")
            for k, v in list(high_miss.items())[:3]:
                lines.append(f"   {k}: {v:.0f}%")
        # 报告路径
        base_dir = os.path.dirname(fpath) if fpath else os.getcwd()
        base_name = os.path.splitext(os.path.basename(fpath))[0] if fpath else "NHANES"
        report_path = os.path.join(base_dir, f"{base_name}_筛选报告.txt")
        if os.path.exists(report_path):
            lines.append(f"\n📄 筛选报告: {base_name}_筛选报告.txt")

        messagebox.showinfo("NHANES 处理完成", "\n".join(lines))

    def _check_db_status_on_startup(self):
        """启动时检查本地数据库是否已存在"""
        db_path = DB_DEFAULT_PATH
        if os.path.exists(db_path) and os.path.getsize(db_path) > 1024:
            self.db_btn.configure(state="disabled", text="✅ 本地数据库已就绪\n(点击重新构建)")
            self.db_query_btn.configure(state="normal")
            self._log("✅ 检测到本地数据库已存在，可直接使用「从本地数据库查询」")
        else:
            self.db_btn.configure(state="normal")
            self.db_query_btn.configure(state="disabled")
            self._log("ℹ️ 本地数据库不存在，如需使用请点击「构建本地数据库」或「导入本地数据库」")

    def _import_db(self):
        """导入已有的本地数据库文件（.db）"""
        p = filedialog.askopenfilename(
            title="选择本地数据库文件",
            filetypes=[("SQLite数据库", "*.db"), ("所有文件", "*.*")],
            initialdir=os.path.dirname(DB_DEFAULT_PATH)
        )
        if not p:
            return
        try:
            import shutil
            shutil.copy2(p, DB_DEFAULT_PATH)
            sz = os.path.getsize(DB_DEFAULT_PATH)
            self.db_btn.configure(state="disabled", text="✅ 本地数据库已就绪\n(点击重新构建)")
            self.db_query_btn.configure(state="normal")
            self._log(f"✅ 数据库已导入: {os.path.basename(p)} ({sz/1024/1024:.1f}MB)")
            messagebox.showinfo("导入成功", f"数据库文件已导入\n大小: {sz/1024/1024:.1f} MB\n可直接使用「从本地数据库查询」")
        except Exception as e:
            messagebox.showerror("导入失败", f"导入数据库失败:\n{e}")
            self._log(f"❌ 导入数据库失败: {e}")

    def _build_db(self):
        """构建本地数据库"""
        if not messagebox.askyesno("确认", "将下载全部13个周期的所有数据\n"
                                          "(共约200-500MB)\n"
                                          "首次构建可能需要较长时间\n"
                                          "确定开始？"):
            return
        self.db_btn.configure(state="disabled")
        self.status_var.set("⏳ 构建本地数据库...")
        self._log("="*50)
        self._log("🗄️ 开始构建本地数据库...")

        self.db_builder = LocalDBBuilder()
        self.db_builder.log_callback = self._log

        thr = threading.Thread(target=self._run_db_build, daemon=True)
        thr.start()
        self.root.after(500, self._check_db_status)

    def _run_db_build(self):
        try:
            result = self.db_builder.build(
                include_mortality=self.mort_var.get(),
                convert_units=self.convert_var.get()
            )
            self.last_db_result = result
        except Exception as e:
            self.last_db_result = {"success": False, "error": str(e)}

    def _check_db_status(self):
        if self.db_builder.is_running:
            self.root.after(1000, self._check_db_status)
        else:
            self.db_btn.configure(state="normal")
            r = self.last_db_result
            if r and r.get("success"):
                self.status_var.set(f"✅ 数据库就绪 ({r['rows']}条)")
                self.db_query_btn.configure(state="normal")
                self._log(f"✅ 本地数据库构建完成: {r['path']}")
                self._log(f"   可随时点击「从本地数据库查询」使用")
            else:
                self.status_var.set("❌ 构建失败")
                self._log(f"❌ 构建失败: {r.get('error', '未知') if r else '无返回'}")

    def _query_db(self):
        """从本地数据库查询"""
        db_path = DB_DEFAULT_PATH
        if not os.path.exists(db_path):
            messagebox.showwarning("提示", "本地数据库不存在，请先构建！")
            return

        win = Tk()
        win.title("从本地数据库查询")
        win.geometry("600x400")

        Label(win, text="NHANES 本地数据库查询", font=("微软雅黑", 12, "bold")).pack(pady=10)
        Label(win, text=f"数据库: {db_path}", font=("微软雅黑", 8), fg="#888").pack()

        Label(win, text="选择周期（可多选，用逗号分隔）\n留空=全部",
              font=("微软雅黑", 9)).pack(pady=5)
        cycle_entry = Entry(win, width=50)
        cycle_entry.insert(0, "2013-2014, 2015-2016")
        cycle_entry.pack(pady=2)

        Label(win, text="选择变量（用逗号分隔）\n留空=全部",
              font=("微软雅黑", 9)).pack(pady=5)
        var_entry = Entry(win, width=50)
        var_entry.insert(0, "性别(Gender), 年龄(Age,岁), BMI(kg/m²), 总胆固醇(TC,mg/dL)")
        var_entry.pack(pady=2)

        def do_query():
            cycles = [c.strip() for c in cycle_entry.get().split(",") if c.strip()]
            vars_list = [v.strip() for v in var_entry.get().split(",") if v.strip()]

            out_path = filedialog.asksaveasfilename(
                defaultextension=".csv",
                filetypes=[("CSV", "*.csv"), ("Excel", "*.xlsx")],
                initialfile="NHANES_query_result.csv"
            )
            if not out_path:
                return

            builder = LocalDBBuilder(db_path)
            builder.log_callback = self._log
            r = builder.query(cycles if cycles else None,
                              vars_list if vars_list else None,
                              out_path)
            if r.get("success"):
                messagebox.showinfo("成功", f"已导出 {r['rows']} 条记录到:\n{out_path}")
                win.destroy()
            else:
                messagebox.showerror("错误", r.get("error", "查询失败"))

        ttk.Button(win, text="查询并导出", command=do_query).pack(pady=10)
        ttk.Button(win, text="取消", command=win.destroy).pack()

    def _preview(self):
        r = self.last_result
        if not r or not r.get("success"):
            messagebox.showinfo("提示", "还没有成功的结果可供预览")
            return
        df = r.get("df")
        if df is None:
            fp = r.get("file_path")
            if not fp or not os.path.exists(fp):
                messagebox.showinfo("提示", "结果文件不存在")
                return
            try:
                df = pd.read_csv(fp, encoding="utf-8-sig") if fp.endswith(".csv") else pd.read_excel(fp)
            except Exception as e:
                messagebox.showerror("错误", f"无法读取文件: {e}")
                return

        self._show_preview_window(df)

    def _show_preview_window(self, df):
        win = Tk()
        win.title("NHANES 数据预览与质控")
        win.geometry("1100x650")

        notebook = ttk.Notebook(win)
        notebook.pack(fill="both", expand=True, padx=5, pady=5)

        # === Tab 1: 数据表格 ===
        tab_data = ttk.Frame(notebook)
        notebook.add(tab_data, text="📊 数据预览")

        toolbar = ttk.Frame(tab_data)
        toolbar.pack(fill="x")
        ttk.Label(toolbar, text=f"共 {len(df)} 行, {len(df.columns)} 列").pack(side="left")

        def _export(fmt):
            name = f"NHANES_export_{datetime.now().strftime('%Y%m%d_%H%M%S')}.{fmt}"
            p = filedialog.asksaveasfilename(defaultextension=f".{fmt}",
                filetypes=[("CSV", "*.csv"), ("Excel", "*.xlsx")], initialfile=name)
            if p:
                (df.to_csv(p, index=False, encoding="utf-8-sig") if fmt == "csv"
                 else df.to_excel(p, index=False, engine="openpyxl"))
                messagebox.showinfo("成功", f"已导出到: {p}")

        ttk.Button(toolbar, text="导出CSV", command=lambda: _export("csv")).pack(side="right", padx=2)
        ttk.Button(toolbar, text="导出Excel", command=lambda: _export("xlsx")).pack(side="right", padx=2)

        frame = ttk.Frame(tab_data)
        frame.pack(fill="both", expand=True)

        tree = ttk.Treeview(frame, show="headings")
        cols = list(df.columns)
        tree["columns"] = cols
        for c in cols:
            tree.heading(c, text=c)
            sample_vals = df[c].dropna().head(3)
            w = max(len(str(c))*10, 80)
            if len(sample_vals) > 0:
                w = max(w, *(len(str(v))*7 for v in sample_vals))
            tree.column(c, width=min(w, 180), anchor="center")

        for _, row in df.head(500).iterrows():
            vals = [str(v) if pd.notna(v) else "" for v in row]
            tree.insert("", END, values=vals)

        hs = ttk.Scrollbar(frame, orient="horizontal", command=tree.xview)
        vs = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        tree.configure(xscrollcommand=hs.set, yscrollcommand=vs.set)
        tree.pack(side="left", fill="both", expand=True)
        vs.pack(side="right", fill="y")
        hs.pack(side="bottom", fill="x")

        # === Tab 2: 质控报告 ===
        tab_qc = ttk.Frame(notebook)
        notebook.add(tab_qc, text="✅ 质控报告")

        qc_r = self.last_result.get("qc_report", {}) if self.last_result else {}
        qc_text = Text(tab_qc, font=("Consolas", 10), wrap="word", padx=10, pady=10)
        qc_text.pack(fill="both", expand=True)

        qc_text.insert(END, "="*60 + "\n")
        qc_text.insert(END, "  NHANES 数据质控报告\n")
        qc_text.insert(END, "="*60 + "\n\n")

        qc_text.insert(END, f"📊 基本信息\n")
        qc_text.insert(END, f"  总记录数: {qc_r.get('总记录数', 'N/A')}\n")
        qc_text.insert(END, f"  总变量数: {qc_r.get('总变量数', 'N/A')}\n\n")

        # 周期分布
        if qc_r.get("周期分布"):
            qc_text.insert(END, f"📅 周期分布\n")
            for k, v in qc_r["周期分布"].items():
                qc_text.insert(END, f"  {k}: {v} 人\n")
            qc_text.insert(END, "\n")

        # 缺失率统计
        if qc_r.get("变量缺失率"):
            qc_text.insert(END, f"⚠️ 变量缺失率(>0%)\n")
            qc_text.insert(END, f"  {'变量名':<30s} {'缺失率':>8s}\n")
            qc_text.insert(END, f"  {'-'*38}\n")
            sorted_miss = sorted(qc_r["变量缺失率"].items(), key=lambda x: -x[1])
            for k, v in sorted_miss:
                bar = "█" * int(v/5) + "░" * (20 - int(v/5))
                qc_text.insert(END, f"  {k:<30s} {v:>6.1f}% {bar}\n")
            qc_text.insert(END, "\n")

        # 数值范围
        if qc_r.get("数值范围"):
            qc_text.insert(END, f"📈 数值变量统计范围\n")
            qc_text.insert(END, f"  {'变量名':<30s} {'最小值':>8s} {'最大值':>8s} {'均值':>8s} {'中位数':>8s}\n")
            qc_text.insert(END, f"  {'-'*62}\n")
            for k, v in qc_r["数值范围"].items():
                qc_text.insert(END, f"  {k:<30s} {v['最小值']:>8.2f} {v['最大值']:>8.2f} {v['均值']:>8.2f} {v['中位数']:>8.2f}\n")

        qc_text.configure(state="disabled")

    # ════════════════════════════════════════════════════════
    # 0.1N 验证回调（严格模式 + 宽松模式）
    # ════════════════════════════════════════════════════════

    def _verify_common(self, strict):
        """0.1N验证通用逻辑"""
        if self.last_result is None or not self.last_result.get("success"):
            messagebox.showwarning("提示", "请先成功导出一次数据！")
            return

        csv_path = self.last_result.get("file_path")
        if not csv_path or not os.path.exists(csv_path):
            messagebox.showwarning("提示", "CSV文件不存在，请重新导出")
            return

        mode_name = "严格" if strict else "宽松"
        if not messagebox.askyesno("确认",
             f"将启动 0.1N {mode_name}验证\n"
             f"  数据: {os.path.basename(csv_path)}\n"
             f"  方法: 标准pandas独立读取CDC原始XPT\n"
             f"  {'容差: 1e-6（需100%精确匹配）' if strict else '容差: 反向换算变量1.0/直接变量1e-6'}\n\n"
             "过程可能需要2-5分钟，确定继续？"):
            return

        self.status_var.set(f"⏳ 0.1N {mode_name}验证中...")
        self._log(f"🔍 启动 0.1N {mode_name}验证...")

        def task():
            try:
                from validation_engine import ValidationEngine
                CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nhanes_cache")
                if not os.path.exists(CACHE_DIR):
                    alt_cache = os.path.join(os.path.expanduser("~"), "nhanes_cache")
                    CACHE_DIR = alt_cache if os.path.exists(alt_cache) else CACHE_DIR

                engine = ValidationEngine(cache_dir=CACHE_DIR)
                result = engine.run(csv_path, sample_ratio=0.1, seed=20260617, strict=strict)
                report_path = engine.export_report()
                self._log(f"✅ {mode_name}验证完成!")
                self._log(f"   通过率: {result['summary']['pass_rate']}")
                self._log(f"   报告: {report_path}")
                messagebox.showinfo(f"{mode_name}验证完成",
                    f"验证方式: 0.1N {mode_name}模式\n"
                    f"总变量: {result['summary']['total_variables_checked']}\n"
                    f"通过: {result['summary']['passed']}\n"
                    f"失败: {result['summary']['failed']}\n"
                    f"通过率: {result['summary']['pass_rate']}\n"
                    f"报告: {os.path.basename(report_path)}")
                self.status_var.set(f"✅ {mode_name}验证: {result['summary']['pass_rate']}")
            except Exception as e:
                self._log(f"❌ {mode_name}验证失败: {e}")
                import traceback
                self._log(traceback.format_exc())
                self.status_var.set("❌ 验证失败")
                messagebox.showerror("错误", str(e))

        threading.Thread(target=task, daemon=True).start()

    def _run_verify_strict(self):
        """严格验证：1e-6 容差"""
        self._verify_common(strict=True)

    def _run_verify_lenient(self):
        """宽松验证：反向换算变量用 0.5 容差"""
        self._verify_common(strict=False)


# ============================================================================
# 主入口
# ============================================================================

def main():
    root = Tk()
    app = NhanesGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
