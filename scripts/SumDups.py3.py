##################################
#                                #
# Last modified 2017/08/07       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string

def run():

    if len(sys.argv) < 4:
        print('usage: python %s input labelfields valuefields outputfilename [-average] [-max]' % sys.argv[0])
        print('       valuefields format: either comma separated, or start:end (including start and end, 0-based)')
        sys.exit(1)
    
    input = sys.argv[1]
    outfilename = sys.argv[4]
    outfile = open(outfilename, 'w')

    fields = sys.argv[2].split(',')
    labelFields=[]
    for f in fields:
        labelFields.append(int(f))
    labelFields.sort()

    print(labelFields)

    doMax = False
    if '-max' in sys.argv:
        doMax = True
        print('will take the max over all instances')

    doAverage = False
    if '-average' in sys.argv:
        doAverage = True
        CountDict={}
        print('will average over instances')

    valueFields=[]
    if ':' in sys.argv[3]:
        fields = sys.argv[3].split(':')
        start = int(fields[0])
        end = int(fields[1])
        for f in range(start,end+1):
            valueFields.append(f)
    else:
        fields = sys.argv[3].split(',')
        for f in fields:
            valueFields.append(int(f))
    valueFields.sort()

    print(valueFields)

    DataDict={}

    LC = 0

    linelist = open(input)
    for line in linelist:
        LC += 1
        if LC % 1000000 == 0:
            print(LC/1000000, 'M lines processed')
        fields = line.replace('\x00','').strip().split('\t')
        if line.startswith('#') or line.startswith('tracking_id'):
            if len(fields) < max(valueFields):
                continue
            outline = '#'
            for ID in labelFields:
                outline = outline + fields[ID] + '\t'
            for ID in valueFields:
                outline = outline + fields[ID] + '\t'
            outfile.write(outline.strip()+'\n')
            continue
        label = []
        for ID in labelFields:
            label.append(fields[ID])
        label = tuple(label)
        if label in DataDict:
            if doAverage:
                CountDict[label]+=1
            for ID in valueFields:
                if DataDict[label][ID] == 'FAIL':
                    continue
                elif fields[ID] == 'FAIL':
                    DataDict[label][ID] = 'FAIL'
                elif doMax:
                    DataDict[label][ID] = max(float(fields[ID]),DataDict[label][ID])
                else:
                    DataDict[label][ID] += float(fields[ID])
        else:
            if doAverage:
                CountDict[label]=1
            DataDict[label]={}
            for ID in valueFields:
                if fields[ID] == 'FAIL':
                    DataDict[label][ID] = 'FAIL'
                else:
                    DataDict[label][ID] = float(fields[ID])

    for label in list(DataDict.keys()):
        label = list(label)
        outline = ''
        for f in label:
            outline = outline + f + '\t'
        label = tuple(label)
        if doAverage:
            for ID in valueFields:
                outline = outline + str(DataDict[label][ID]/CountDict[label]) + '\t'
        else:
            for ID in valueFields:
                outline = outline + str(DataDict[label][ID]) + '\t'
        outfile.write(outline.strip()+'\n') 
            
    outfile.close()
   
run()

##################################
#                                #
# Last modified 2017/08/07       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string

def run():

    if len(sys.argv) < 4:
        print('usage: python %s input labelfields valuefields outputfilename [-average] [-max]' % sys.argv[0])
        print('       valuefields format: either comma separated, or start:end (including start and end, 0-based)')
        sys.exit(1)
    
    input = sys.argv[1]
    outfilename = sys.argv[4]
    outfile = open(outfilename, 'w')

    fields = sys.argv[2].split(',')
    labelFields=[]
    for f in fields:
        labelFields.append(int(f))
    labelFields.sort()

    print(labelFields)

    doMax = False
    if '-max' in sys.argv:
        doMax = True
        print('will take the max over all instances')

    doAverage = False
    if '-average' in sys.argv:
        doAverage = True
        CountDict={}
        print('will average over instances')

    valueFields=[]
    if ':' in sys.argv[3]:
        fields = sys.argv[3].split(':')
        start = int(fields[0])
        end = int(fields[1])
        for f in range(start,end+1):
            valueFields.append(f)
    else:
        fields = sys.argv[3].split(',')
        for f in fields:
            valueFields.append(int(f))
    valueFields.sort()

    print(valueFields)

    DataDict={}

    LC = 0

    linelist = open(input)
    for line in linelist:
        LC += 1
        if LC % 1000000 == 0:
            print(LC/1000000, 'M lines processed')
        fields = line.replace('\x00','').strip().split('\t')
        if line.startswith('#') or line.startswith('tracking_id'):
            if len(fields) < max(valueFields):
                continue
            outline = '#'
            for ID in labelFields:
                outline = outline + fields[ID] + '\t'
            for ID in valueFields:
                outline = outline + fields[ID] + '\t'
            outfile.write(outline.strip()+'\n')
            continue
        label = []
        for ID in labelFields:
            label.append(fields[ID])
        label = tuple(label)
        if label in DataDict:
            if doAverage:
                CountDict[label]+=1
            for ID in valueFields:
                if DataDict[label][ID] == 'FAIL':
                    continue
                elif fields[ID] == 'FAIL':
                    DataDict[label][ID] = 'FAIL'
                elif doMax:
                    DataDict[label][ID] = max(float(fields[ID]),DataDict[label][ID])
                else:
                    DataDict[label][ID] += float(fields[ID])
        else:
            if doAverage:
                CountDict[label]=1
            DataDict[label]={}
            for ID in valueFields:
                if fields[ID] == 'FAIL':
                    DataDict[label][ID] = 'FAIL'
                else:
                    DataDict[label][ID] = float(fields[ID])

    for label in list(DataDict.keys()):
        label = list(label)
        outline = ''
        for f in label:
            outline = outline + f + '\t'
        label = tuple(label)
        if doAverage:
            for ID in valueFields:
                outline = outline + str(DataDict[label][ID]/CountDict[label]) + '\t'
        else:
            for ID in valueFields:
                outline = outline + str(DataDict[label][ID]) + '\t'
        outfile.write(outline.strip()+'\n') 
            
    outfile.close()
   
run()

