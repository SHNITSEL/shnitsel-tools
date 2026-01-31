_datatree_level_attribute_key = "DataTree_Level"

# Map to get the appropriate corresponding class name for each
# level name in a hierarchical tree structure.
# Has some legacy names mapped to new type names
DataTreeLevelMap = {
    "root": "ShnitselDBRoot",
    "compound": "CompoundGroup",
    "group": "DataGroup",
    "trajectory": "DataLeaf",
    "data": "DataLeaf",
    "ShnitselDBRoot": "ShnitselDBRoot",
    "CompoundGroup": "CompoundGroup",
    "TrajectoryGroup": "DataGroup",
    "TrajectoryData": "DataLeaf",
    "DataGroup": "DataGroup",
    "DataLeaf": "DataLeaf",
}
