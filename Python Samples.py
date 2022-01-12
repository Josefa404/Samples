#Beispiel 1: Funktion um mathematische Funktionen zu konkatenieren und der relativen Fehler aus der Berechnung zu bestimmen

import math
import numpy as np

print("2. a)")

#concat-Funktion:
def concat(g1, g2):
    def comp(a):
        return g2(g1(a))
    return comp

#Test

def f1(n):
    return 2*n

def f2(n):
    return n + 3

f = concat(f1, f2)

print(f(42))

print("2. c)")

#Funktion um schnell den relativen Fehler bestimmen zu können:
def relativeError(x,y):
    return np.abs(x-y)/np.abs(x)

#f(x)
def ff(n):
    return math.sqrt(n+1)-1

#f^-1(x)
def ffhme(n):
    return n*(n+2)

def generateData(x):
    v1 = concat(ff, ffhme)
    val1 = v1(x)
    v2 = concat(ffhme, ff)
    val2 = v2(x)
    err1 = relativeError(x,val1)
    err2 = relativeError(x,val2)
    return (val1,val2,err1,err2)

#Beispiel 2: Berechnung der relativen Kondition von Elementaroperationen:

def subtract(x1, x2):
    # berechnet die relative Konditionszahl der Funktion x1 - x2
    val = x1 - x2
    cond = np.abs(x1) + (np.abs(x2)) / np.abs(x1-x2)
    return val, cond

def times2(x):
    # berechnet die relative Konditionszahl der Funktion 2*x
    val = 2*x
    cond = 2
    return val, cond

#Beispiel 3: Implementierung von Auswertungsbäumen 

def visualize_tree(b):
    for pre, fill, node in RenderTree(b):
        if callable(node.name):
            printname = node.name.__name__
        else: printname = node.name
        print("%s%s" % (pre, printname))

bF2 = Node(add) #f5
be = Node(subtract, parent = bF2) #f1
bf = Node(square, parent = bF2) #f2
b04 = Node(x2, parent = bf) #x2
ba = Node(square, parent = be) #f1
b01 = Node(x1, parent = ba)  #x1
bb = Node(times2, parent = be) #f4
bc = Node(multiply, parent = bb) #f3
b02 = Node(x1, parent = bc) #x1
b03 = Node(x2, parent = bc) #x2

print("Visualisierung vom Auswertungsbaum:")
visualize_tree(bF2)

