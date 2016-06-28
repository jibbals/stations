# -*- coding: utf-8 -*-

from datetime import datetime, timedelta

def tau_to_date(tau):
    '''
        Convert tau list into datetime list
    ''' 
    date0 = datetime(1985,1,1,0,0)
    newtau = tau.astype(float)
    dates= [ date0+timedelta(hours=h) for h in newtau ]
    return dates
