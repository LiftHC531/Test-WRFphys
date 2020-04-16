#!/bin/bash
find -name "*.F" | xargs grep -iR "WDM6" > find.txt
