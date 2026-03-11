#!/bin/bash
# V3 外部工具安装脚本
# 安装 Merizo（结构域分割）、HMMER（Pfam 搜索）、Pfam-A 数据库
#
# 使用方法:
#   conda activate cryo-pipeline
#   bash scripts/install_tools.sh
#
# 需要网络代理时先执行:
#   source <(curl -sSL http://deploy.i.h.pjlab.org.cn/infra/scripts/setup_proxy.sh)

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
TOOLS_DIR="$PROJECT_ROOT/tools"

mkdir -p "$TOOLS_DIR"

echo "=============================="
echo "V3 外部工具安装"
echo "=============================="

# 1. PyTorch CPU-only（Merizo 依赖）
echo ""
echo "[1/4] 安装 PyTorch CPU-only..."
if python -c "import torch" 2>/dev/null; then
    echo "  PyTorch 已安装，跳过"
else
    pip install torch --index-url https://download.pytorch.org/whl/cpu || {
        echo "  [WARN] PyTorch 安装失败，Merizo 将不可用"
    }
fi

# 2. Merizo（结构域分割）
echo ""
echo "[2/4] 安装 Merizo..."
MERIZO_DIR="$TOOLS_DIR/merizo"
if [ -d "$MERIZO_DIR" ] && [ -f "$MERIZO_DIR/predict.py" ]; then
    echo "  Merizo 已存在: $MERIZO_DIR"
else
    git clone https://github.com/psipred/merizo.git "$MERIZO_DIR" || {
        echo "  [WARN] Merizo 克隆失败，domain 标签将全零"
    }
fi

# Merizo 依赖
if [ -d "$MERIZO_DIR" ]; then
    pip install einops networkx "rotary-embedding-torch>=0.3.0" natsort 2>/dev/null || {
        echo "  [WARN] Merizo 依赖安装部分失败"
    }
fi

# 3. HMMER（Pfam 搜索）
echo ""
echo "[3/4] 安装 HMMER..."
if command -v hmmscan &>/dev/null; then
    echo "  HMMER 已安装: $(which hmmscan)"
else
    conda install -c bioconda hmmer -y 2>/dev/null || {
        echo "  [WARN] HMMER conda 安装失败，尝试手动下载..."
        # 回退：下载静态二进制
        HMMER_URL="http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz"
        HMMER_TMP="$TOOLS_DIR/hmmer_build"
        mkdir -p "$HMMER_TMP"
        wget -q -O "$HMMER_TMP/hmmer.tar.gz" "$HMMER_URL" && \
            tar xzf "$HMMER_TMP/hmmer.tar.gz" -C "$HMMER_TMP" --strip-components=1 && \
            cd "$HMMER_TMP" && ./configure --prefix="$TOOLS_DIR/hmmer_install" && \
            make -j4 && make install && \
            export PATH="$TOOLS_DIR/hmmer_install/bin:$PATH" || {
                echo "  [WARN] HMMER 手动安装也失败，Pfam split 将回退到 MMseqs2"
            }
        cd "$PROJECT_ROOT"
    }
fi

# 4. Pfam-A 数据库
echo ""
echo "[4/4] 下载 Pfam-A 数据库..."
PFAM_DIR="$TOOLS_DIR/pfam"
mkdir -p "$PFAM_DIR"
if [ -f "$PFAM_DIR/Pfam-A.hmm.h3m" ]; then
    echo "  Pfam-A 数据库已就绪"
else
    PFAM_URL="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    if [ ! -f "$PFAM_DIR/Pfam-A.hmm" ]; then
        echo "  下载 Pfam-A.hmm.gz (~2GB)..."
        wget -q --show-progress -P "$PFAM_DIR" "$PFAM_URL" || {
            echo "  [WARN] Pfam-A 下载失败（需要代理？），Pfam split 将回退到 MMseqs2"
        }
        if [ -f "$PFAM_DIR/Pfam-A.hmm.gz" ]; then
            echo "  解压中..."
            gunzip "$PFAM_DIR/Pfam-A.hmm.gz"
        fi
    fi
    if [ -f "$PFAM_DIR/Pfam-A.hmm" ]; then
        echo "  压缩索引 (hmmpress)..."
        hmmpress "$PFAM_DIR/Pfam-A.hmm" || {
            echo "  [WARN] hmmpress 失败"
        }
    fi
fi

echo ""
echo "=============================="
echo "安装完成！验证状态："
echo "=============================="
echo -n "  PyTorch: "; python -c "import torch; print(torch.__version__)" 2>/dev/null || echo "未安装"
echo -n "  Merizo:  "; [ -f "$TOOLS_DIR/merizo/predict.py" ] && echo "OK ($TOOLS_DIR/merizo)" || echo "未安装"
echo -n "  HMMER:   "; command -v hmmscan 2>/dev/null || echo "未安装"
echo -n "  Pfam-A:  "; [ -f "$PFAM_DIR/Pfam-A.hmm.h3m" ] && echo "OK ($PFAM_DIR)" || echo "未就绪"
echo ""
