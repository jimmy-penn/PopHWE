import re
import os
import mimetypes


#
# Get time difference
#
def GetTimeDifference(t1, t2):
    dsec = 0
    dusec = 0
    dhour = 0
    dmin = 0

    fdt = t2 - t1

    dusec = int((fdt - int(fdt)) * 1000000)
    dt = int(fdt)
    
    if dt >= 3600:
        dhour = int(dt/3600.0)
        dt -= dhour * 3600
        
    if dt >= 60:
        dmin = int(dt/60.0)
        dt -= dmin * 60
        
    dsec = dt

    return dhour, dmin, dsec, dusec        


#
# Get time difference in string
#
def GetTimeDiffString(t1, t2):
    dhour, dmin, dsec, dusec = GetTimeDifference(t1, t2)
    
    dstr = f''
    if dhour > 0:
        dstr += f' {dhour} hour'
        
        if dhour > 1:
            dstr += f's'
          
    if dmin > 0:
        dstr += f' {dmin} minute'
        
        if dmin > 1:
            dstr += f's'

    if dsec > 0:
        dstr += f' {dsec} second'
        
        if dsec > 1:
            dstr += f's'

    if dhour == 0 and dmin == 0 and dusec > 0:
        dstr += f' {dusec} microsecond'
        
        if dusec > 1:
            dstr += f's'

    dstr = dstr.strip()

    return dstr
