package com.aliyun.odps.exttable.bam;

import java.util.*;

public class Deduper {
    private Set<Integer> cacheSet;

    public Deduper(int maxSize)
    {
        this.cacheSet = Collections.newSetFromMap(new LinkedHashMap<Integer, Boolean>(maxSize * 4 / 3 + 1)
        {
            @Override
            protected boolean removeEldestEntry(final Map.Entry eldest)
            {
                return size() > maxSize;
            }
        });
    }

    public boolean isDeplicated(int val) {
        return this.cacheSet.contains(val);
    }

    public void put(int val)
    {
        this.cacheSet.add(val);
    }
}
