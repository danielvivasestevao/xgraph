# Xgraph.py
The core of this project.
It contains three classes: Xgraph, XgraphList and XgraphGUI.
Importing the file with `from Graphs.Xgraph import *` imports Xgraph and
XgraphList.

## Xgraph objects
A wrapper for `networkx.Graph` objects which allows calculation of coexistence and redundancy properties and graph clustering on the graph object. Can be saved to and loaded from binary files (using `pickle`) with built-in methods.

### Parameters
`Xgraph(
gps, lognormal, path_loss_exponent, variance,
path_loss_one_meter, poisson_point, transmit_power,
lam, dim_x, dim_y)`

### Saving and loading an Xgraph
Save an Xgraph object to a binary file called xgraph. Instead of just the file
name, you can use the whole path to the file. (You should __always__ use the
absolute path to a file when saving it. See The troubleshooting section below for
more information.)

`xg.save("/path/to/file/xgraph")`

You can load an Xgraph file with the static class method `Xgraph.load()`.

`xg = Xgraph.load("/path/to/file/xgraph")` 

### Calculating coexistence and redundancy
Calculate coexistence and redundancy of an Xgraph.

If you pass the functions an integer parameter k, coexistence/redundancy will be calculated up to a path length of k hops; otherwise, it is calculated up to the graph's diameter.

`# xg is an Xgraph object
xg.calc_coexistence()
xg.calc_redundancy()
xg.calc_coexistence(k=5)
xg.calc_redundancy(k=2)`

Calling these functions on an Xgraph object will calculate a `Coexistence` or a `Redundancy` object respectively and assign it to the Xgraph. You can either get the exact values in form of a Python dictionary using
`xg.get_coexistence()._dict` or `xg.get_redundancy()._dict`, or you can get the statistics of the respective object in form of a Python dictionary with `xg.get_coexistence_stats(percent)` or `xg.get_redundancy_stats(percent)`, where `percent` is `True` for percentage values or `False` for absolute values.

Verbose example:
```
>>> xg = Xgraph(lognormal=True, path_loss_exponent=3.5,
                 variance=10, path_loss_one_meter=38.0,
                 poisson_point=True, transmit_power=-10.0,
                 lam=100, dim_x=100)
>>> xg.calc_coexistence()
>>> xg.get_coexistence_stats(percent=True)
{'strong': {1: 0.9214758751182592, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0, 6: 1.0, 7: 1.0}, 'weak': {1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0, 6: 1.0, 7: 1.0}}
>>> xg.get_coexistence_stats(percent=False)
{'strong': {1: 974, 2: 1057, 3: 1057, 4: 1057, 5: 1057, 6: 1057, 7: 1057}, 'weak': {1: 1057, 2: 1057, 3: 1057, 4: 1057, 5: 1057, 6: 1057, 7: 1057}}
```

### Clustering

#### Cluster types
(At the moment, only one naive clustering method is implemented. More methods will be added in the future.)

`naive`: Picks a random node of the graph, declares it as a cluster head, and adds all nodes in its `k`-hop neighbourhood to its cluster. Repeat until all nodes are either cluster heads or belong to a cluster. If and only if two nodes of different clusters are connected, their cluster heads are connected.

#### Usage
You can cluster an Xgraph object using the `cluster(cluster_type: str, k: int)` method. `cluster_type` is the name of the clustering method, and `k` is the range of the cluster. For example, for `cluster_type="naive"` and `k=3`, each cluster would consist of nodes connected by at most three hops.

`# xg is an Xgraph object
xg.cluster("naive", 3)`

The clustered graph is added to the Xgraph's `_cluster` dictionary. An Xgraph can have one clustered version of itself for every pair of `cluster_type` and `k` value. You can access a clustered graph with the `get_cluster(cluster_type: str, k: int=None)` method.

`clustered_graph = xg.get_cluster("naive", k)`

If you cluster an Xgraph more than once using the same parameters, the old clustered graph will be overwritten.

Clustered graphs are Xgraph objects themselves, which means you can calculate coexistence and redundancy on them, save them as seperate files, cluster them further, etc. 

## XgraphList objects
XgraphList is used to store multiple Xgraph objects. It is an extension of the
Python `list` type. It only accepts Xgraph objects as elements and supports
functions for calculating various average values for all contained Xgraph
objects. It is also possible to calculate redundancy and coexistence
dictionaries of all graphs in a list in multiprocessing mode.

### Adding Xgraph objects to an XgraphList
Since XgraphList inherits Python's `list` type, it behaves like a `list` object (for the most part). You can add an Xgraph object by appending it to the list as follows.

`# xg is an Xgraph object, xgl is an XgraphList object
xgl.append(xg)`

In __"small memory mode"__ (see the chapter below), you can also append an Xgraph's ID as a `str` or `UUID` object to an XgraphList. If the object you are trying to append to the list is not of a supported type, an `Exception` will be raised.

### Saving and loading an XgraphList
Analogous to saving and loading an Xgraph file. In 

### Calculating coexistence and redundancy
It is possible to calculate coexistence/redundancy of all Xgraph objects in an XgraphList in one go. The following functions exploit multicore processors in order to calculate coexistence/redundancy of multiple Xgraph objects at once.

If you pass the functions an integer parameter k, coexistence/redundancy will be calculated up to a path length of k hops; otherwise, coexistence/redundancy up to the diameter of the Xgraph object with the biggest diameter out of all Xgraph objects in the list is calculated.

`# xgl is an XgraphList object
xgl.calc_coexistence()
xgl.calc_redundancy()
xgl.calc_coexistence(k=5)
xgl.calc_redundancy(k=2)`

### Calculating average/maximum values
Functions like `avrg_diameter()`, `avrg_node_degree()`, `avrg_node_number()` or `max_diameter()` call the respective function on all Xgraph objects in the list and calculate the average/maximum of the returned values.

### "Small memory" mode
If an XgraphList object is initialized with `xgl = XgraphList(small_mem=True)`, Xgraph objects appent to it will not be loaded into memory, but saved as a binary file on disk. Instead of the whole Xgraph object, just its ID will be appent to the XgraphList as a reference to the binary file, and the Xgraph object itself is only loaded into memory whenever a computation on it has to be done.

Note: When in small memory mode, the XgraphList and its Xgraph objects will be saved in the current working directory. To delete them, use `xgl._delete_from_disk()`. This will delete the XgraphList file as well as the directory where all of its Xgraph files are written to. **Please consider that the whole directory will be deleted**, so do not move any other files into it. The directorie's name is "Xgraph_files_<ID of the XgraphList>".


## Troubleshooting

### Import errors
If you are getting import errors when trying to execute a script, you might
have to add the project's root folder (in this case that would be `redcoe`) to
your Python path.
Just execute `export PYTHONPATH=${PYTHONPATH}:/path/to/project/folder` to do so.
To add it permanently to your Python path variable, add that line to your
`~/.bashrc`.

### Xgraph(List) objects are not saved
Try saving the object using the __absolute__ path instead of a relative one.