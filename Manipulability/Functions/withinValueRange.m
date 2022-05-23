function bool = withinValueRange(value, comparisonValue, range)
    belowRange = value <= comparisonValue+range;
    aboveRange = value >= comparisonValue-range;
    bool = belowRange && aboveRange;
end