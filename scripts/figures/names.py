"""One display name per network, used by every figure, table and data file of the paper.

A network keeps its model id where that is a BiGG-style id. Three model files are named otherwise and get a short name:
Chlorella_variabilis_iAJ526 and Saha2011iRS1597 their BiGG-style ids iAJ526 and iRS1597, ArabidopsisCoreModel the
name AraCore. The three life stages of iIS312 are abbreviated the same way everywhere (Try, Ama, Epi for the
trypomastigote, amastigote and epimastigote models), because the figures have no room for the full names.

    from names import network      # network("iIS312_Amastigote") == "iIS312 Ama"
"""
NAME = {"ArabidopsisCoreModel": "AraCore",
        "Chlorella_variabilis_iAJ526": "iAJ526",
        "Saha2011iRS1597": "iRS1597",
        "iIS312_Trypomastigote": "iIS312 Try",
        "iIS312_Amastigote": "iIS312 Ama",
        "iIS312_Epimastigote": "iIS312 Epi"}


def network(m):
    """display name of the model file id m"""
    return NAME.get(m, m)
