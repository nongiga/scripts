function str=myExtractBefore(str, endStr)
    toRem=contains(str, endStr);
    str(toRem)=extractBefore(str(toRem),endStr);
end