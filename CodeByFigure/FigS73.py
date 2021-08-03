#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 12:12:32 2021

@author: ryan
"""
import numpy as np
import matplotlib.pyplot as plt
from random import random
from functools import reduce
import tensorflow as tf
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)







def translist(list_here):
    return [list(x) for x in zip(*list_here)]

def flatten(t):
    return [item for sublist in t for item in sublist]

def Check(x):
    return False == all(map(lambda y: all(map(lambda z: z == 0, y)), x))

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109503
genesRaw = [i for i in open('ChondrogenesisOnlineData.txt',"r")]

donor = [0,1,2]
Time = [0,1,3,7,14,21]
TestDonor = 0


NotTestDonor = slice(0+(TestDonor==0),3-(TestDonor==2),1+(TestDonor==1))

cats = genesRaw[0].split()
geneind = cats.index('geneName')
gene = [i.split()[geneind] for i in genesRaw[1:]]
timeind = translist([[cats.index(i)-2 for i in cats if 'D'+str(j)+'D' in i] for j in [1,2,3]])

exp = [i.split()[2:] for i in genesRaw[1:]]
exp2 = tuple(map(lambda  x: [(float(x[k])) for k in translist(timeind)[TestDonor]], exp))
exp = tuple(map(lambda  x: [[(float(x[k])) for k in l] for l in translist(timeind)[NotTestDonor]], exp))


exp = [[gene[i],exp[i],exp2[i]] for i in range(len(exp))]
exp = list(filter(lambda x: Check(x[1]), exp))
# exp.sort(key = lambda x: max( [(max(i)-min(i))/(min(i)) for i in translist(x[1])] ), reverse = False)
# exp = exp[:1000]

gene = tuple(map(lambda x: x[0], exp)) 
exp2 = tuple(map(lambda  x: x[2], exp)) 
exp = tuple(map(lambda  x: x[1], exp)) 

Normize = lambda x,y: np.array(x)/np.array(y)
exp = [Normize(i,exp[gene.index('EEF2')]) for i in exp]
exp2 = [Normize(i,exp2[gene.index('EEF2')]) for i in exp2]


Exp = {gene[i]: exp[i] for i in range(len(gene))}
TestData = {gene[i]: exp2[i] for i in range(len(gene))} 
AllExp = {gene[i]: exp[i] + [exp2[i]] for i in range(len(gene))} 

e = 2.7182818284590452353602874713527 
def getProbFunc(geneSampled,ValSeen):
    LikelyhoodNotT = tuple(map(lambda x: abs(sum(map(lambda y: y - ValSeen, x))), translist(Exp[geneSampled])))
    
    ProbT = tuple(map(lambda x: 0.01 + reduce(lambda x,y: x*y, [LikelyhoodNotT[i] for i in range(len(Time)) if (i != x)]), range(len(Time))))
    Porportionality = 1/sum(ProbT)
    return tuple(map(lambda x: x*Porportionality, ProbT))

VarianceFromNoise = 10
def getProbFuncBayes(geneSampled,ValSeen):
    GeneData = translist(Exp[geneSampled])
    # GeneData = translist(AllExp[geneSampled])
    mu = tuple(map(lambda x: sum(x)/len(x), GeneData))
    numtop = tuple(map(lambda x, m: sum(map(lambda y: (y - m)**2,x)), GeneData, mu))
    sig = tuple(map(lambda x, m: (x/(len(numtop)-1) + VarianceFromNoise)**0.5, numtop, mu))
    ProbT = tuple(map(lambda m, s: 0.01 + e**(-0.5*((m - ValSeen)/s)**2)/s, mu, sig))
    Porportionality = 1/sum(ProbT)
    return tuple(map(lambda x: x*Porportionality, ProbT))


def TestSpecifiedGenes(MeasuredAtNow,geneSampled):
    AllProbs = tuple(map(lambda x: getProbFunc(x, TestData[x][MeasuredAtNow]), geneSampled))
    probs = reduce(lambda x,y: tuple(map(lambda a, b: a*b/sum(x) , x,y)), AllProbs)
    probs = tuple(map(lambda x: x/sum(probs), probs))
    error = sum(map(lambda x, y: x*y, probs, Time)) - Time[MeasuredAtNow]
    return error

def TestSpecifiedGenesBayes(MeasuredAtNow,geneSampled):
    AllProbs = tuple(map(lambda x: getProbFuncBayes(x, TestData[x][MeasuredAtNow]), geneSampled))
    probs = reduce(lambda x,y: tuple(map(lambda a, b: a*b/sum(x) , x,y)), AllProbs)
    probs = tuple(map(lambda x: x/sum(probs), probs))
    error = sum(map(lambda x, y: x*y, probs, Time)) - Time[MeasuredAtNow]
    return error


def TestSpecifiedGenesArgmax(MeasuredAtNow,geneSampled):
    AllProbs = tuple(map(lambda x: getProbFunc(x, TestData[x][MeasuredAtNow]), geneSampled))
    probs = reduce(lambda x,y: tuple(map(lambda a, b: a*b/sum(x) , x,y)), AllProbs)
    probs = tuple(map(lambda x: x/sum(probs), probs))
    error = Time[np.argmax(probs)] - Time[MeasuredAtNow]
    return error

def TestSpecifiedGenesBayesArgmax(MeasuredAtNow,geneSampled):
    AllProbs = tuple(map(lambda x: getProbFuncBayes(x, TestData[x][MeasuredAtNow]), geneSampled))
    probs = reduce(lambda x,y: tuple(map(lambda a, b: a*b/sum(x) , x,y)), AllProbs)
    probs = tuple(map(lambda x: x/sum(probs), probs))
    error = Time[np.argmax(probs)] - Time[MeasuredAtNow]
    return error

"""
0:11 for 05 genes, double check 50 times, 4000 datas*num_genes
1:16 for 20 genes, double check 50 times, 4000 datas*num_genes 0:20 for 1000
7:00 for 50 genes, double check 50 times, 4000 datas*num_genes 1:37 for 1000
"""
# num_genes = 12
for num_genes in (5,12,20,50):
    double_check_number = 50
    error_set = []
    error_setb = []
    for iter_double_check in range(double_check_number):
        GenesToDeduceTime = random.sample(gene, num_genes)
        # GenesToDeduceTime = ['ACTB', 'COL1A1', 'EEF2', 'GAPDH', 'RUNX1', 'SOX9', 'SPP1']
        TrainingInfo = translist([translist(Exp[i]) for i in GenesToDeduceTime])
        TrainingDay = np.array(Time)
        TestInfo = np.array(translist([TestData[i] for i in GenesToDeduceTime]))
        TestDay = np.array(Time)
        
        
        def make1data():
            ThisTimeInd = int(random.random()*len(Time))
            thisRand = random.random()
            # Input = [np.random.normal(j[0] + thisRand*(j[1]-j[0]), np.std(j)*(1-abs(thisRand-0.5)*2), 1) for j in TrainingInfo[ThisTimeInd]]
            # Input = [np.random.normal(j[0] + thisRand*(j[1]-j[0]), np.std(j)*(abs(thisRand-0.5)*2), 1) for j in TrainingInfo[ThisTimeInd]]
            Input = [np.random.normal(j[0] + thisRand*(j[1]-j[0]), np.std(j), 1) for j in TrainingInfo[ThisTimeInd]]
            Output = ThisTimeInd
            return [Input, Output]
        
        
        train_data = translist([make1data() for i in range(1000*num_genes)])
        x_train = np.array(train_data[0])
        y_train = np.array(train_data[1])
        
        
        test_output = [[t] for t in range(len(Time))]
        test_input = [[[i] for i in TestInfo[t]] for t in range(len(Time))]
        x_test = np.array(test_input)
        y_test = np.array(test_output)
        
        model = tf.keras.models.Sequential([
          tf.keras.layers.Flatten(input_shape= x_train[0].shape ),
          tf.keras.layers.Dense(20, activation='relu'),
          tf.keras.layers.Dropout(0.2),
          tf.keras.layers.Dense(len(Time), activation='softmax')
        ])
        predictions = model(x_train[:1]).numpy()
        predictions
        tf.nn.softmax(predictions).numpy()
        loss_fn = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)
        loss_fn(y_train[:1], predictions).numpy()
        
        model.compile(optimizer='adam',
                      loss=loss_fn,
                      metrics=['accuracy'])
        model.fit(x_train, y_train,validation_split = 0.1, epochs=20)
        model.evaluate(x_test,  y_test, verbose=2)
        
        probability_model = tf.keras.Sequential([
          model,
          tf.keras.layers.Softmax()
        ])
        probability_model(x_test)
        
        
        MyGeneErrorsNeuralNet = tuple(np.dot(model.predict(x_test),Time) - Time)
        MyGeneErrors = tuple([TestSpecifiedGenes(t,GenesToDeduceTime) for t in range(len(Time))])
        MyGeneErrorsBayes = tuple([TestSpecifiedGenesBayes(t,GenesToDeduceTime) for t in range(len(Time))])
        
        MyGeneErrorsNeuralNetArgmax = tuple([Time[i]-Time[j] for j,i in enumerate(np.argmax(model.predict(x_test), axis=1))])
        MyGeneErrorsArgmax = tuple([TestSpecifiedGenesArgmax(t,GenesToDeduceTime) for t in range(len(Time))])
        MyGeneErrorsBayesArgmax = tuple([TestSpecifiedGenesBayesArgmax(t,GenesToDeduceTime) for t in range(len(Time))])
        
        if iter_double_check == 0: #Check for bias in estimator
            e1b = [i/double_check_number for i in MyGeneErrorsNeuralNet]
            e2b = [i/double_check_number for i in MyGeneErrors]
            e3b = [i/double_check_number for i in MyGeneErrorsBayes]  
            e4b = [i/double_check_number for i in MyGeneErrorsNeuralNetArgmax]
            e5b = [i/double_check_number for i in MyGeneErrorsArgmax]
            e6b = [i/double_check_number for i in MyGeneErrorsBayesArgmax]            
        else:
            e1b = [i/double_check_number+j for i,j in zip(MyGeneErrorsNeuralNet,e1b)]
            e2b = [i/double_check_number+j for i,j in zip(MyGeneErrors,e1b)]
            e3b = [i/double_check_number+j for i,j in zip(MyGeneErrorsBayes,e1b)]
            e4b = [i/double_check_number+j for i,j in zip(MyGeneErrorsNeuralNetArgmax,e1b)]
            e5b = [i/double_check_number+j for i,j in zip(MyGeneErrorsArgmax,e1b)]
            e6b = [i/double_check_number+j for i,j in zip(MyGeneErrorsBayesArgmax,e1b)]
        
        if iter_double_check == 0: #Check standard deviation
            e1 = [abs(i)/double_check_number for i in MyGeneErrorsNeuralNet]
            e2 = [abs(i)/double_check_number for i in MyGeneErrors]
            e3 = [abs(i)/double_check_number for i in MyGeneErrorsBayes]  
            e4 = [abs(i)/double_check_number for i in MyGeneErrorsNeuralNetArgmax]
            e5 = [abs(i)/double_check_number for i in MyGeneErrorsArgmax]
            e6 = [abs(i)/double_check_number for i in MyGeneErrorsBayesArgmax]            
        else:
            e1 = [abs(i)/double_check_number+j for i,j in zip(MyGeneErrorsNeuralNet,e1)]
            e2 = [abs(i)/double_check_number+j for i,j in zip(MyGeneErrors,e1)]
            e3 = [abs(i)/double_check_number+j for i,j in zip(MyGeneErrorsBayes,e1)]
            e4 = [abs(i)/double_check_number+j for i,j in zip(MyGeneErrorsNeuralNetArgmax,e1)]
            e5 = [abs(i)/double_check_number+j for i,j in zip(MyGeneErrorsArgmax,e1)]
            e6 = [abs(i)/double_check_number+j for i,j in zip(MyGeneErrorsBayesArgmax,e1)]
        print('Runthrough: '+str(iter_double_check+1))
    error_setb.extend([e1b, e2b, e3b, e4b, e5b, e6b])
    error_set.extend([e1, e2, e3, e4, e5, e6])
    
    
    
    fig, ax = plt.subplots(figsize=(12, 9))
    im = ax.imshow(error_set, cmap = 'RdBu')
    ax.set_xticks(np.arange(len(error_set[0])))
    ax.set_yticks(np.arange(len(error_set)))
    ax.set_xticklabels(['Day '+str(i) for i in Time])
    ax.set_yticklabels(['Neural Network Probability','MLP Probability','Bayes Probability','Neural Network Argmax','MLP Argmax','Bayes Argmax'])
    ax.set_title("Error magnatude of differentiation state predictor: "+str(num_genes)+' genes')
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('|Error| [Days]', rotation=-90, va="bottom")
    plt.savefig('EEF2MultiplePredictors'+str(num_genes)+'genesAbsval.png')
    
    
    
    fig, ax = plt.subplots(figsize=(12, 9))
    im = ax.imshow(error_setb, cmap = 'RdBu')
    ax.set_xticks(np.arange(len(error_setb[0])))
    ax.set_yticks(np.arange(len(error_setb)))
    ax.set_xticklabels(['Day '+str(i) for i in Time])
    ax.set_yticklabels(['Neural Network Probability','MLP Probability','Bayes Probability','Neural Network Argmax','MLP Argmax','Bayes Argmax'])
    ax.set_title("Error of differentiation state predictor: "+str(num_genes)+' genes')
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('Error [Days]', rotation=-90, va="bottom")
    plt.savefig('EEF2MultiplePredictors'+str(num_genes)+'genes.png')


