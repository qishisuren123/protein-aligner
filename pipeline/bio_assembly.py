"""
模块: Bio Assembly - 生物学组装体展开

功能：
- 将不对称单元展开为完整的生物学组装体
- 例如 apoferritin（24聚体）从 ~3000 原子展开到 ~76000 原子
- 处理无组装体信息的回退（返回原始结构）
- 安全限制防止内存溢出
"""
import logging
import gemmi

logger = logging.getLogger(__name__)


def expand_to_assembly(structure, assembly_id="1", max_chains=200):
    """
    将结构展开为完整的生物学组装体

    参数:
        structure: gemmi.Structure 对象
        assembly_id: 组装体 ID（默认 "1"，即第一个组装体）
        max_chains: 最大链数限制，防止内存溢出

    返回:
        展开后的 gemmi.Structure（如果无组装体信息则返回原始结构的副本）
    """
    # 检查是否有组装体信息
    if not structure.assemblies:
        logger.info("  无组装体信息，使用原始结构")
        return structure.clone()

    # 查找指定的组装体
    assembly = None
    for asm in structure.assemblies:
        if asm.name == assembly_id:
            assembly = asm
            break

    if assembly is None:
        # 尝试使用第一个组装体
        assembly = structure.assemblies[0]
        logger.info(f"  未找到组装体 '{assembly_id}'，使用第一个: '{assembly.name}'")

    # 预估展开后的链数
    n_orig_chains = len(structure[0])
    n_ops = sum(len(gen.operators) for gen in assembly.generators)
    estimated_chains = n_ops * n_orig_chains if n_ops > 0 else n_orig_chains
    logger.info(f"  组装体 '{assembly.name}': {n_ops} 个对称操作, "
                f"预估 {estimated_chains} 条链")

    if estimated_chains > max_chains:
        logger.warning(f"  展开后链数 ({estimated_chains}) 超过限制 ({max_chains})，"
                       f"使用原始结构")
        return structure.clone()

    # 执行展开
    try:
        model = structure[0]
        new_st = structure.clone()
        new_model = gemmi.make_assembly(
            assembly, model, gemmi.HowToNameCopiedChain.Short
        )
        new_st[0] = new_model

        # 统计展开结果
        n_new_chains = len(new_st[0])
        n_new_atoms = sum(
            1 for chain in new_st[0]
            for res in chain if res.name != "HOH"
            for _ in res
        )
        logger.info(f"  展开完成: {n_orig_chains} -> {n_new_chains} 条链, "
                     f"{n_new_atoms} 个原子")
        return new_st
    except Exception as e:
        logger.error(f"  组装体展开失败: {e}")
        return structure.clone()


def get_assembly_info(structure):
    """
    获取结构的组装体信息摘要

    返回:
        dict: 组装体信息
    """
    info = {
        "n_assemblies": len(structure.assemblies),
        "assemblies": [],
    }
    for asm in structure.assemblies:
        n_ops = sum(len(gen.operators) for gen in asm.generators)
        info["assemblies"].append({
            "name": asm.name,
            "n_operators": n_ops,
        })
    return info
