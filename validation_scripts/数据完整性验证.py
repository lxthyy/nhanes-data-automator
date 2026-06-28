# -*- coding: utf-8 -*-
"""
NHANES 涓嬭浇鍣?鏁版嵁瀹屾暣鎬ч獙璇佽剼鏈?浣跨敤鏂规硶: python 鏁版嵁瀹屾暣鎬ч獙璇?py
杈撳嚭: 楠岃瘉缁撴灉鎵撳嵃鍒版帶鍒跺彴

楠岃瘉鍘熺悊锛氳鍙栦笅杞藉櫒杈撳嚭鐨凜SV锛屾寜璋冩煡鍛ㄦ湡鍒嗙粍锛?閫愬垪璁＄畻缂哄け鐜囷紝鏍囪"鍛ㄦ湡鏃犳暟鎹?锛堣鍛ㄦ湡100%缂哄け锛変笌"涓綋缂哄け"銆?"""
import pandas as pd
import sys

sys.stdout.reconfigure(encoding='utf-8')
CSV = r'..\NHANES_E_F_G_瀹屾暣鏁版嵁_v205.csv'
df = pd.read_csv(CSV, encoding='utf-8-sig', low_memory=False)

print('鏁版嵁瀹屾暣鎬ч獙璇佺粨鏋?)
print(f'鎬绘牱鏈噺: {len(df)} 浜?)
print(f'鎬诲彉閲忔暟: {len(df.columns)} 鍒?)
print(f'鍛ㄦ湡: {df["璋冩煡鍛ㄦ湡"].value_counts().to_dict()}')
print()
print('=' * 70)

for cycle in df['璋冩煡鍛ㄦ湡'].unique():
    subset = df[df['璋冩煡鍛ㄦ湡'] == cycle]
    print(f'\n--- {cycle} (N={len(subset)}) ---')

    total_vars = 0
    high_na = []  # >70% 缂哄け
    mid_na = []   # 30-70%
    zero_na = 0   # 鍛ㄦ湡缂哄け锛堟暣鍒楀叏NA锛?
    for col in df.columns:
        if col in ('搴忓彿(SEQN)', '璋冩煡鍛ㄦ湡'):
            continue
        total_vars += 1
        na_pct = subset[col].isna().mean() * 100

        if na_pct == 100:
            zero_na += 1
        elif na_pct > 70:
            high_na.append((col, na_pct))
        elif na_pct > 30:
            mid_na.append((col, na_pct))

    print(f'  鎬诲彉閲? {total_vars}')
    print(f'  鉁?浣庣己澶?<=30%): '
          f'{total_vars - len(high_na) - len(mid_na) - zero_na}')
    print(f'  馃煛 涓害缂哄け(30-70%): {len(mid_na)}')
    print(f'  馃敶 楂樼己澶?>70%): {len(high_na)}')

    if zero_na > 0:
        print(f'  鉂?鍛ㄦ湡鏃犳妫€娴?100%缂哄け): {zero_na}')

    if high_na:
        print(f'  楂樼己澶卞彉閲忎妇渚?')
        for col, pct in sorted(high_na, key=lambda x: -x[1])[:5]:
            n_valid = subset[col].notna().sum()
            print(f'    {col}: {pct:.1f}% (鏈夋晥N={n_valid})')

print()
print('娉? 缂哄け鐜?100%鐨勫垪琛ㄧず璇ュ彉閲忓湪璇ュ懆鏈熸棤妫€娴嬮」鐩紝鏄疌DC璁捐濡傛锛岄潪鏁版嵁閿欒銆?)
input('\n鎸?Enter 閿€€鍑?..')
