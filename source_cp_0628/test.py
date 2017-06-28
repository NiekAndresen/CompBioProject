import re
import numpy as np
example = 'N(0.05)L(0.05)G(0.05)Ksda-loop(0.05)V(0.05)G(0.05)S(0.29)K(0.28)Ccm(0.05)Ccm(0.05)K'
def get_neighborhood_list(linkMap):
    matches = re.findall(r'(\w+(?:-\w+)?)\((\d\.\d{2})\)', linkMap)
    aaList = []; valList = []; posList = []; origPosList = np.arange(len(matches)); valSum = 0
    for i,match in enumerate(matches):
        value = float(match[1])
        if len(match[0])==1:
            aaList += [match[0]]
            valList += [value]
            posList += [origPosList[i]]
        elif match[0][-1].isupper():
            aaList += [match[0][-1]]
            valList += [value]
            posList += [origPosList[i]]
        valSum += value
    valList = np.array(valList) / valSum
    print(posList)
    posList -= posList[np.argmax(valList)]
    return aaList, valList, posList

match = get_neighborhood_list(example)
print(match)


#next: see how the entries in the linkmap conincide with positions in the neighborhood


