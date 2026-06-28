# -*- coding: utf-8 -*-
"""
NHANES 涓嬭浇鍣?鍒嗙被鍙橀噺鏄犲皠楠岃瘉鑴氭湰
浣跨敤鏂规硶: python 鍒嗙被鍙橀噺鏄犲皠楠岃瘉.py
杈撳嚭: 楠岃瘉缁撴灉鎵撳嵃鍒版帶鍒跺彴

楠岃瘉鍘熺悊锛氳鍙栦笅杞藉櫒杈撳嚭鐨凜SV锛屽鍒嗙被鍙橀噺鍒楁鏌ュ敮涓€鍊硷紝
纭鎵€鏈夊€煎潎涓轰腑鏂?鑻辨枃鏍囩锛屾棤CDC鍘熷鏁板瓧缂栫爜娈嬬暀銆?"""
import pandas as pd
import sys

sys.stdout.reconfigure(encoding='utf-8')
CSV = r'..\NHANES_E_F_G_瀹屾暣鏁版嵁_v205.csv'

df = pd.read_csv(CSV, encoding='utf-8-sig', low_memory=False)

# 寰呴獙璇佺殑鍏抽敭鍒嗙被鍙橀噺
key_cats = [
    '鎬у埆(Gender)', '绉嶆棌(Race)', '鏁欒偛绋嬪害(Education)', '濠氬Щ鐘跺喌(Marital)',
    '绯栧翱鐥呰瘖鏂?DM-Dx)', '鏇惧惛100鏀儫(Smok100+)',
    '鏇鹃ギ閰?2娆′互涓?EverDrink)', '鑷瘎鍋ュ悍(Health-Gen)',
    '鍝枠(Asthma)', '鎵撻季(Snore)',
    '鏇捐瘖鏂珮琛€鍘?HTN-Dx)', '鏈嶇敤闄嶅帇鑽?HTN-Meds)',
    '鐫＄湢闅滅(Sleep-Dx)', '鐫＄湢闅滅(Sleep-Prob)',
    '鍙備笌鐘舵€?Status)', '鎬€瀛曠姸鎬?Pregnant)',
    '浣跨敤澶勬柟鑽?Rx-Use)', '鍑虹敓鍥?Birth-US)', '鍏皯(Citizen)',
]

pass_n = 0
fail_n = 0
print('鍒嗙被鍙橀噺鏄犲皠楠岃瘉缁撴灉:')
print('=' * 60)

for label in key_cats:
    found_cols = [c for c in df.columns if label in c]
    if not found_cols:
        print(f'  鈿狅笍 鏈壘鍒? {label}')
        continue
    col = found_cols[0]
    vals = sorted(df[col].dropna().unique())
    # 妫€鏌ユ槸鍚︽湁鏁板瓧缂栫爜娈嬬暀
    has_digits = all(
        str(v).replace('.', '').replace('-', '').replace(' ', '').isdigit()
        for v in vals if pd.notna(v)
    )
    if has_digits:
        status = '鉂?鍚暟瀛楃紪鐮?
        fail_n += 1
    else:
        status = '鉁?
        pass_n += 1
    print(f'  {status} {label}: {vals}')

print()
print(f'缁撴灉: {pass_n} 閫氳繃, {fail_n} 澶辫触, '
      f'閫氳繃鐜?{pass_n/(pass_n+fail_n)*100:.0f}%' if pass_n + fail_n > 0 else '鏃犲彉閲忓彲姣斿')
input('\n鎸?Enter 閿€€鍑?..')
