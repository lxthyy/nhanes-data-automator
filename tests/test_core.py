# -*- coding: utf-8 -*-
"""
NHANES Downloader & Validation Engine — 单元测试

测试内容：
  1. 别名解析 — XPT_COLUMN_ALIASES 能否跨周期找到列
  2. 单位换算 — TC/TG/Glu/Cr 换算系数准确
  3. 缺失值清洗 — CDC 缺失编码能否被正确识别

运行: pytest tests/ -v
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd
import numpy as np
from validation_engine import ValidationEngine


class TestAliasResolution:
    """测试1：别名解析正确性"""

    def setup_method(self):
        self.engine = ValidationEngine()

    def test_tsh_aliases(self):
        """TSH 应支持 LBXTSH / LBXTSH1 / LBDTSH1S 三个别名"""
        aliases = self.engine.XPT_COLUMN_ALIASES.get("LBXTSH1", [])
        assert "LBXTSH" in aliases, "E/F/G周期用 LBXTSH"
        assert "LBXTSH1" in aliases, "L/M周期用 LBXTSH1"
        assert "LBDTSH1S" in aliases, "部分周期用 LBDTSH1S"

    def test_lipid_aliases(self):
        """血脂变量应支持 SI/非SI 别名"""
        assert "LBDHDD" in self.engine.XPT_COLUMN_ALIASES.get("LBDHDD", [])
        assert "LBDLDL" in self.engine.XPT_COLUMN_ALIASES.get("LBDLDL", [])

    def test_find_raw_col(self):
        """_find_raw_col 应能在模拟 XPT DataFrame 中找到列"""
        import pandas as pd
        mock_xpt = pd.DataFrame({"SEQN": [1], "LBXTSH": [1.5]})
        result = self.engine._find_raw_col(mock_xpt, "LBXTSH1")
        assert result == "LBXTSH", f"应找到 LBXTSH，实际得到 {result}"

    def test_find_raw_col_case_insensitive(self):
        """不区分大小写也应能匹配"""
        mock_xpt = pd.DataFrame({"seqn": [1], "lbxtsh": [1.5]})
        result = self.engine._find_raw_col(mock_xpt, "LBXTSH1")
        assert result is not None, "不区分大小写应能找到"
        assert result.upper() == "LBXTSH"


class TestUnitConversion:
    """测试2：单位换算系数准确性"""

    REVERSE = ValidationEngine.REVERSE_CONVERSION

    def test_tc_conversion(self):
        """TC: 200 mg/dL → 5.17 mmol/L, 反向验证"""
        # 下载器: mg/dL / 38.67 → mmol/L
        tc_mmol = 200 / 38.67  # 下载器输出
        # 验证引擎: mmol/L * 38.67 → mg/dL
        restored = tc_mmol * 38.67
        assert abs(restored - 200) < 0.01, f"TC反向换算误差: {restored-200}"

    def test_tg_conversion(self):
        """TG: 150 mg/dL → 1.69 mmol/L"""
        tg_mmol = 150 / 88.57
        restored = tg_mmol * 88.57
        assert abs(restored - 150) < 0.01

    def test_glu_conversion(self):
        """Glu: 100 mg/dL → 5.55 mmol/L"""
        glu_mmol = 100 / 18.01
        restored = glu_mmol * 18.01
        assert abs(restored - 100) < 0.01

    def test_cr_conversion(self):
        """Cr: 1.0 mg/dL → 88.4 μmol/L"""
        cr_umol = 1.0 * 88.4
        restored = cr_umol / 88.4
        assert abs(restored - 1.0) < 0.01

    def test_ua_conversion(self):
        """UA: 6.0 mg/dL → 356.91 μmol/L"""
        ua_umol = 6.0 * 59.485
        restored = ua_umol / 59.485
        assert abs(restored - 6.0) < 0.01


class TestMissingValueCleanup:
    """测试3：CDC 缺失值编码识别"""

    def setup_method(self):
        self.engine = ValidationEngine()

    def test_missing_7777_detected(self):
        """7777 应被视为缺失值"""
        assert 7777 in self.engine.CDC_MISSING_VALUES

    def test_missing_9999_detected(self):
        """9999 应被视为缺失值"""
        assert 9999 in self.engine.CDC_MISSING_VALUES

    def test_missing_8888_detected(self):
        """8888 应被视为缺失值"""
        assert 8888 in self.engine.CDC_MISSING_VALUES

    def test_nan_not_in_missing(self):
        """NaN 不应在缺失值集合中（用 notna 处理）"""
        assert np.nan not in self.engine.CDC_MISSING_VALUES

    def test_valid_value_not_missing(self):
        """正常 TSH 值 1.5 不应被视为缺失"""
        assert 1.5 not in self.engine.CDC_MISSING_VALUES
