function ret = getCurPrecise(curGen,maxGen)

initVal = 2;
endVal = 6;

ret = round(initVal + (endVal- initVal) * curGen / maxGen);