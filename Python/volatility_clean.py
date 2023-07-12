
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.interpolate
import glob
import random
import matplotlib
import utility

task='ana' #'ana' for anagrams;  'cra' for compound remote association
excludeID,excludeTrials=utility.getExcludes(task)
includeID=utility.getIncludeID(task)
prefix='log_'



outputFolder = '/projects/p31274/Drexel/processed/'+task+'_chan/'

sourceFolder = '/projects/p31274/Drexel/export/' + task + '_chan_allFreq/'
allchannels = [x.lower() for x in utility.channellist]

chanindex = list(range(1,len(allchannels)+1)) #use all channels

if task=='ana':
    timeToStimOffset=.5 #anagram onset is .5 sec too early
else:
    timeToStimOffset=0
includePreOnset=-1 #sec
def parseFilename(filename):
    value=filename.split('.')[0]
    values=value.split('_')
    subid=values[1]
    type=values[2]
    index=values[3]
    return subid,type,int(index)
def calcVolatility(subid,prefix,sourceFolder,includeComps):
    volDataForSub=pd.DataFrame()

    for file in glob.glob(sourceFolder +subid + '/'+prefix+subid+'*'):
        filename=file.split('/')[-1]
        subid, type, index = parseFilename(filename)
        if (subid,type,index) in excludeTrials:continue
        allData,featLbl=utility.preprocess(file,includeComps,includeFreqs=utility.includeFreqs_full)
        allData['timeToStim']=allData['timeToStim']-timeToStimOffset
        duration=allData.iloc[0]['timeToStim']-allData.iloc[0]['timeToResp']
        nDataPoints = len(
            allData[(allData['timeToStim'] < 0) & (allData['timeToStim'] > includePreOnset)][featLbl].diff())
        volDF0=allData[(allData['timeToStim']<0)&(allData['timeToStim']>includePreOnset)][featLbl].diff().std()

        volDF0_co = allData[(allData['timeToStim'] < -.25) & (allData['timeToStim'] > includePreOnset-.25)][featLbl].diff().std()


        volDF1=allData[(allData['timeToStim']>0)&(allData['timeToStim']<duration/2)][featLbl].diff().std()


        volDF2 = allData[(allData['timeToResp'] < 0) & (allData['timeToStim'] > duration / 2)][featLbl].diff().std()

        volDF_l2s = allData[(allData['timeToResp'] < 0) & (allData['timeToResp'] > -2)][featLbl].diff().std()


        labelUnpack = [x.split(',') for x in featLbl]
        channellables=[utility.channellist[int(np.float(x[0]))-1] for x in labelUnpack]

        freqlabel = [x[1] for x in labelUnpack]

        volDF=pd.DataFrame({'pre':volDF0,'precutoff': volDF0_co,'early':volDF1,'late':volDF2,'last2s':volDF_l2s,
                            'chan':channellables,'freq':freqlabel,'npts':nDataPoints})
        volDF['subid']=subid
        volDF['type']=type
        volDF['index']=index

        volDataForSub = pd.concat([volDataForSub,volDF])

    return volDataForSub

allVols=pd.DataFrame()
alltrials=pd.DataFrame()
for subid in includeID:
    trialinfo=pd.read_csv(sourceFolder +subid + '/trials.csv',index_col=None)
    trialinfo['subid']=subid
    alltrials = pd.concat([alltrials,trialinfo])
    if subid in excludeID:continue
    volDataForSub = calcVolatility(subid,prefix,sourceFolder,chanindex)
    allVols=pd.concat([allVols,volDataForSub])

allVols.to_csv(outputFolder+'trialsummary_volatility_allchan.csv',index=None)
