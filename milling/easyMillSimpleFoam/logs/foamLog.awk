# Awk script for OpenFOAM log file extraction
BEGIN {
    Iteration=0
    resetCounters()
}

# Reset counters used for variable postfix
function resetCounters() {
    SeparatorCnt=0
    TimeCnt=0
    UxCnt=0
    UxFinalResCnt=0
    UxItersCnt=0
    UyCnt=0
    UyFinalResCnt=0
    UyItersCnt=0
    UzCnt=0
    UzFinalResCnt=0
    UzItersCnt=0
    clockTimeCnt=0
    contCumulativeCnt=0
    contGlobalCnt=0
    contLocalCnt=0
    epsAvgCnt=0
    epsMaxCnt=0
    epsMinCnt=0
    epsilonCnt=0
    epsilonFinalResCnt=0
    epsilonItersCnt=0
    executionTimeCnt=0
    kCnt=0
    kAvgCnt=0
    kFinalResCnt=0
    kItersCnt=0
    kMaxCnt=0
    kMinCnt=0
    pCnt=0
    pFinalResCnt=0
    pItersCnt=0
    # Reset counters for 'Solving for ...'
    for (varName in subIter)
    {
        subIter[varName]=0
    }
}

# Extract value after columnSel
function extract(inLine,columnSel,outVar,a,b)
{
    a=index(inLine, columnSel)
    b=length(columnSel)
    split(substr(inLine, a+b),outVar)
    gsub("[,:]","",outVar[1])
}

# Iteration separator (increments 'Iteration')
/^[ \t]*Time = / {
    Iteration++
    resetCounters()
}

# Time extraction (sets 'Time')
/^[ \t]*Time = / {
    extract($0, "Time = ", val)
    Time=val[1]
}

# Skip whole line with singularity variable
/solution singularity/ {
    next;
}

# Extract: 'Solving for ...'
/Solving for/ {
    extract($0, "Solving for ", varNameVal)

    varName=varNameVal[1]
    file=varName "_" subIter[varName]++
    file="logs/" file
    extract($0, "Initial residual = ", val)
    print Time "\t" val[1] > file

    varName=varNameVal[1] "FinalRes"
    file=varName "_" subIter[varName]++
    file="logs/" file
    extract($0, "Final residual = ", val)
    print Time "\t" val[1] > file

    varName=varNameVal[1] "Iters"
    file=varName "_" subIter[varName]++
    file="logs/" file
    extract($0, "No Iterations ", val)
    print Time "\t" val[1] > file
}

# Extract: 'Separator'
/^[ \t]*Time = / {
    extract($0, "Time = ", val)
    file="logs/Separator_" SeparatorCnt
    print Time "\t" val[1] > file
    SeparatorCnt++
}

# Extract: 'Time'
/^[ \t]*Time = / {
    extract($0, "Time = ", val)
    file="logs/Time_" TimeCnt
    print Time "\t" val[1] > file
    TimeCnt++
}

# Extract: 'clockTime'
/ClockTime = / {
    extract($0, "ClockTime =", val)
    file="logs/clockTime_" clockTimeCnt
    print Time "\t" val[1] > file
    clockTimeCnt++
}

# Extract: 'contCumulative'
/time step continuity errors :/ {
    extract($0, "cumulative = ", val)
    file="logs/contCumulative_" contCumulativeCnt
    print Time "\t" val[1] > file
    contCumulativeCnt++
}

# Extract: 'contGlobal'
/time step continuity errors :/ {
    extract($0, " global = ", val)
    file="logs/contGlobal_" contGlobalCnt
    print Time "\t" val[1] > file
    contGlobalCnt++
}

# Extract: 'contLocal'
/time step continuity errors :/ {
    extract($0, "sum local = ", val)
    file="logs/contLocal_" contLocalCnt
    print Time "\t" val[1] > file
    contLocalCnt++
}

# Extract: 'epsAvg'
/bounding epsilon,/ {
    extract($0, "average: ", val)
    file="logs/epsAvg_" epsAvgCnt
    print Time "\t" val[1] > file
    epsAvgCnt++
}

# Extract: 'epsMax'
/bounding epsilon,/ {
    extract($0, "max: ", val)
    file="logs/epsMax_" epsMaxCnt
    print Time "\t" val[1] > file
    epsMaxCnt++
}

# Extract: 'epsMin'
/bounding epsilon,/ {
    extract($0, "min: ", val)
    file="logs/epsMin_" epsMinCnt
    print Time "\t" val[1] > file
    epsMinCnt++
}

# Extract: 'executionTime'
/ExecutionTime = / {
    extract($0, "ExecutionTime = ", val)
    file="logs/executionTime_" executionTimeCnt
    print Time "\t" val[1] > file
    executionTimeCnt++
}

# Extract: 'kAvg'
/bounding k,/ {
    extract($0, "average: ", val)
    file="logs/kAvg_" kAvgCnt
    print Time "\t" val[1] > file
    kAvgCnt++
}

# Extract: 'kMax'
/bounding k,/ {
    extract($0, "max: ", val)
    file="logs/kMax_" kMaxCnt
    print Time "\t" val[1] > file
    kMaxCnt++
}

# Extract: 'kMin'
/bounding k,/ {
    extract($0, "min: ", val)
    file="logs/kMin_" kMinCnt
    print Time "\t" val[1] > file
    kMinCnt++
}

# End
