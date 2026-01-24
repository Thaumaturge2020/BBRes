# Code and Datasets of BBRes for Maximum-k-Defective Clique Searching.

Quick start:

```
cd BBRes
make
./BBRes -g graph_path -k k_num
```

## A. Input graph data format:

### Format: */edges.txt

text graph:

+ first line: $n$,$m$

+ then:$m$ lines,each line contains two integers $u$,$v$ as point label.

Our point labels are range from 0 to n-1.

## B. Algorithm BBRes

the whole procedure for searching Maximum k-Defective Clique is located at [BBRes/](https://github.com/Thaumaturge2020/BBRes/tree/main/BBRes)

### 1. Compile

```
make
```

### 2. Run

```
./BBRes -g graph_path -k k_num
```

### 3.An example

```
cd BBRes
make
./BBRes -g ../data/Ca-GrQc -k 3
```

### 4.About major components in codes

+ BBRes corresponds to (./BBRes/kDefective_BB_bitset.hpp,./BBRes/kDefective_BB_matrix.hpp);
+ PreProcessing corresponds to (./BBRes/Graph.cpp)
# BBRes
