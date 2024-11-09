#<<<<< Elastic neutral axis >>>>>
import math
from math import pi
from math import ceil
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, r2_score
from matplotlib.colors import LinearSegmentedColormap
from statistics import mean
from sympy import symbols, solve
from sympy.abc import b,k
from collections import OrderedDict
import pandas as pd
import csv
#tempdist=[[900,0],[180,0.5],[90,1]]
#temp=[400,500,600,700,800,900]
temp=[500]
lamda_group=[]
lamda_amb_group=[]
M_b_fi_group=[]
tempdist_group=[]
chi_LT_fi_group=[]
chi_LT_group=[]
length=[11000,8000,6000,3500]
# length=[6000]
length_amb=[17000,16500,16000,15500,15000,14500,14000,13500,13000,12500,12000,11500,11000,10500,10000,9500,9000,8500,8000,7500,7000,6500,6000,5500,5000,4500,4000,3500,3000,2500,2000,1500]
def y_neutral_stre(y_neutral,stre_dist):
    for i in range(len(stre_dist)-1):
        if stre_dist[i][1]<=y_neutral and stre_dist[i+1][1]>y_neutral:
            y_neutral_stre=(y_neutral-stre_dist[i][1])*(stre_dist[i+1][0]-stre_dist[i][0])/(stre_dist[i+1][1]-stre_dist[i][1])+stre_dist[i][0]
    return y_neutral_stre
def interpolation(a,stre_dist):
    for i in  range(len(stre_dist)-1):
        if stre_dist[i][1] <= a and stre_dist[i+1][1] > a:
            y_inter = (a-stre_dist[i][1])*(stre_dist[i+1][0]-stre_dist[i][0])/(stre_dist[i+1][1]-stre_dist[i][1])+stre_dist[i][0]
    return y_inter
def ela_red(temp):
    coefficient=[[20,1],[100,1],[200,0.9],[300,0.8],[400,0.7],[500,0.6],[600,0.31],[700,0.13],[800,0.09],[900,0.0675],[1000,0.045],[1100,0.0225],[1200,0]]
    for i in range(len(coefficient)-1):
        if temp >= coefficient[i][0] and temp < coefficient[i+1][0]:
            co_of_red = (temp-coefficient[i][0])/(coefficient[i+1][0]-coefficient[i][0])*(coefficient[i+1][1]-coefficient[i][1])+coefficient[i][1]
    return co_of_red
def intemp(a,b,c):
    d_temp=(a-b[0])*(b[1]-c[1])/(b[0]-c[0])+b[1]
    return d_temp
def Sort(sub_li):
    sub_li.sort(key = lambda x: x[1])
    return sub_li

def elastic_neutral_axis(d,w,tf,tw,tempdist,fpl,E,Iz,Iw,It,Wpl,G):
    temp=[]
    int_temp=[]
    rev_int_temp=[]
    # temperatures
    for i in range(len(tempdist)):
        temp.append(tempdist[i][0])
    temp.sort()
    for i in range(len(temp)-1):
        for j in range(1,13):
            if 100*j<tempdist[i][0] and 100*j>tempdist[i+1][0] or 100*j>tempdist[i][0] and 100*j<tempdist[i+1][0]:
                int_temp.append(100*j)
    rev_int_temp=int_temp.copy()
    rev_int_temp.reverse()
    # Interpolation of temperatures and positions
    for k in range(len(tempdist)+len(int_temp)-1):
        if tempdist[k][0]>tempdist[k+1][0]:
            for j in range(1,13):
                if tempdist[k][0]>j*100 and tempdist[k+1][0]<j*100:
                    pos=intemp(j*100,tempdist[k],tempdist[k+1])
                    tempdist.insert(k+1,[j*100,pos])
        elif tempdist[k][0]<tempdist[k+1][0]:
            for j in range(1,13):
                if tempdist[k][0]<j*100 and tempdist[k+1][0]>j*100:
                    pos=intemp(j*100,tempdist[k],tempdist[k+1])
                    tempdist.insert(k+1,[j*100,pos])
    print('temp_dist=',tempdist)
    disp_tempdist=[]
    for i in tempdist:
        disp_tempdist.append([i[0],i[1]])
    ela_dist=[]
    i=0
    for i in range(len(tempdist)):
        ela_dist.append(tempdist[i].copy())
    i=0
    for ela_pt in ela_dist:
        co_of_red=ela_red(ela_pt[0])
        ela_dist[i][0]=co_of_red
        i=i+1
        
    #add=[]
    #for i in range(len(tempdist)-1):
    #    num1=(abs(tempdist[0][0]-tempdist[-1][0])//100)
    #    num=int(abs(tempdist[0][0]-tempdist[-1][0])//100)
    #    if num1 >1:
    #        if tempdist[-1][0]>tempdist[0][0]:
    #            for j in range(1,num+1):
    #                temp=tempdist[0][0]+j*100
    #                pos=intemp(temp,tempdist[i],tempdist[i+1])
    #                add.append([temp,pos])
    #        if tempdist[-1][0]<tempdist[0][0]:
    #            for j in range(1,num+1):
    #                temp=tempdist[0][0]-j*100
    #                pos=intemp(temp,tempdist[i],tempdist[i+1])
    #                add.append([temp,pos])  
    #temp_set=[]
    #for i in add:
    #   for j in tempdist:
    #       temp_set.append(j[0])
    #   if i[0] not in temp_set:
    #       tempdist.append(i)

    
    tempdist=Sort(tempdist)
    print(tempdist)
    #print('<<<<<<  Plastic Bending Moment Calculation  >>>>>>')
    #print('\n<<<<<<  Section Geometry  >>>>>>')
    #print('Depth  d = '+str(d)+' mm')
    #print('Width  w = '+str(w)+' mm')
    #print('Flange thickness  tf = '+str(tf)+' mm')
    #print('Web thickness  tw = '+str(tw)+' mm')
    #print('\n<<<<<<  Section Properties  >>>>>>')
    #print('Second moment of inertia of weak axis Iz = '+str(round(Iz*1e-4,2))+' cm^4')
    #print('Warping constant Iw = '+str(round(Iw*1e-6,2))+' cm^6')
    #print('Torsional constant It = '+str(round(It*1e-4,2))+' cm^4')
    #print('Plastic section modulus Wpl = '+str(round(Wpl*1e-3,2))+' cm^3')
    #print('\n<<<<<<  Material properties  >>>>>>')
    #print('Yielding strength f_y = '+str(fpl)+' N/mm^2')
    #print('Elastic modulus E = '+str(E)+' N/mm^2')
    #print('Shear modulus G = '+str(G)+' N/mm^2')

    #print('\n<<<<<<  Temperature Distribution  >>>>>>')
    #print('Temp.distribution = [[Temperature,Position (x/d, x is the distance from top of the beam)],...]')
    disp_tempdist=[]
    for i in tempdist:
        disp_tempdist.append([i[0],i[1]])
    #print('Temp.distribution =',disp_tempdist)
    ela_dist=[]
    i=0
    for i in range(len(tempdist)):
        ela_dist.append(tempdist[i].copy())
    i=0
    # Strength Reduction according to EC3
    i=0
    for ela_pt in ela_dist:
        co_of_red=ela_red(ela_pt[0])
        ela_dist[i][0]=co_of_red
        i=i+1
    a_s=tf/d
    b_s=(d-tf)/d
    uf_lower_ela=interpolation(a_s,ela_dist)
    lf_upper_ela=interpolation(b_s,ela_dist)
    uf_lower_temp=interpolation(a_s,tempdist)
    lf_upper_temp=interpolation(b_s,tempdist)    
    num_d_uf=0
    num_d_lf=0
    num_temp_uf=0
    num_temp_lf=0
    for i in range(len(tempdist)):
        if tempdist[i][1]==a_s:
            num_temp_uf+=1
    for i in range(len(tempdist)):
        if tempdist[i][1]==b_s:
            num_temp_lf+=1
    if num_temp_uf==0:
        tempdist.append([uf_lower_temp,a_s])
    if num_temp_lf==0:
        tempdist.append([lf_upper_temp,b_s])
    tempdist=Sort(tempdist)  
    for i in range(len(ela_dist)):
        if ela_dist[i][1]==a_s:
            num_d_uf+=1 
    for i in range(len(ela_dist)):
        if ela_dist[i][1]==b_s:
            num_d_lf+=1
    if num_d_uf==0:
        ela_dist.append([uf_lower_ela,a_s])
    if num_d_lf==0:
        ela_dist.append([lf_upper_ela,b_s])
    ela_dist=Sort(ela_dist)
    # ela_dist.append([uf_lower_ela,a_s])
    # ela_dist.append([lf_upper_ela,b_s])
    #for i in range(len(ela_dist)-1):
    #    k=0
    #    for j in range(len(ela_dist)-1):
    #        if ela_dist[i]==ela_dist[j]:
    #            k+=1
    #    if k==1:
    #        ela_dist_new.append(ela_dist[i])
    #ela_dist=ela_dist_new
    
    print('tempdist=',tempdist)
    print('ela_dist=',ela_dist)
    counter_half=0
    area1=0
    area2=0
    area3=0
    area4=0
    y_neu=[]
    y_neutral=0
    for i in range(len(ela_dist)-1):
        if ela_dist[i+1][1]==0.5:
            break
        else:
            counter_half=counter_half+1
    i=0
    disp_ela_dist=[]
    for i in ela_dist:
        disp_ela_dist.append([i[0],i[1]])
    
    for i in range(len(ela_dist)): 
        area1=0
        area2=0
        area3=0
        area4=0 
        area_t=0
        area_b=0    
        for j in range(i):
            if ela_dist[j+1][1]<=tf/d :
                area1=area1+((ela_dist[j][0]-(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*ela_dist[i][1]*(ela_dist[j+1][1]-ela_dist[j][1])+0.5*(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)*ela_dist[i][1]-0.5*(ela_dist[j][0]-(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)-1/3*(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**3-ela_dist[j][1]**3))*w
                #area1=area1+(ela_dist[j][0]*((ela_dist[j][1]-ela_dist[i][1])**2)**0.5+ela_dist[j+1][0]*((ela_dist[j+1][1]-ela_dist[i][1])**2)**0.5)*(ela_dist[j+1][1]-ela_dist[j][1])/2*w
            elif ela_dist[j][1]>=(d-tf)/d and ela_dist[j+1][1]<=1:
                area1=area1+((ela_dist[j][0]-(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*ela_dist[i][1]*(ela_dist[j+1][1]-ela_dist[j][1])+0.5*(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)*ela_dist[i][1]-0.5*(ela_dist[j][0]-(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)-1/3*(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**3-ela_dist[j][1]**3))*w
            elif ela_dist[j][1]>=tf/d and ela_dist[j+1][1]<=(d-tf)/d:
                area1=area1+((ela_dist[j][0]-(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*ela_dist[i][1]*(ela_dist[j+1][1]-ela_dist[j][1])+0.5*(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)*ela_dist[i][1]-0.5*(ela_dist[j][0]-(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)-1/3*(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**3-ela_dist[j][1]**3))*tw
                #area1=area1+(ela_dist[j][0]*((ela_dist[j][1]-ela_dist[i][1])**2)**0.5+ela_dist[j+1][0]*((ela_dist[j+1][1]-ela_dist[i][1])**2)**0.5)*(ela_dist[j+1][1]-ela_dist[j][1])/2*tw 
        for k in range(i,len(ela_dist)-1):
            if ela_dist[k+1][1]<=tf/d :
                area2=area2+(0.5*(ela_dist[k][0]-(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)+1/3*(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**3-ela_dist[k][1]**3)-(ela_dist[k][0]-(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*ela_dist[i][1]*(ela_dist[k+1][1]-ela_dist[k][1])-0.5*(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)*ela_dist[i][1])*w
                #area2=area2+(ela_dist[k][0]*((ela_dist[k][1]-ela_dist[i][1])**2)**0.5+ela_dist[k+1][0]*((ela_dist[k+1][1]-ela_dist[i][1])**2)**0.5)*(ela_dist[k+1][1]-ela_dist[k][1])/2*w
            elif ela_dist[k][1]>=(d-tf)/d and ela_dist[k+1][1]<=1:
                area2=area2+(0.5*(ela_dist[k][0]-(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)+1/3*(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**3-ela_dist[k][1]**3)-(ela_dist[k][0]-(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*ela_dist[i][1]*(ela_dist[k+1][1]-ela_dist[k][1])-0.5*(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)*ela_dist[i][1])*w
            elif ela_dist[k][1]>=tf/d and ela_dist[k+1][1]<=(d-tf)/d:
                area2=area2+(0.5*(ela_dist[k][0]-(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)+1/3*(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**3-ela_dist[k][1]**3)-(ela_dist[k][0]-(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*ela_dist[i][1]*(ela_dist[k+1][1]-ela_dist[k][1])-0.5*(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)*ela_dist[i][1])*tw
                #area2=area2+(ela_dist[k][0]*((ela_dist[k][1]-ela_dist[i][1])**2)**0.5+ela_dist[k+1][0]*((ela_dist[k+1][1]-ela_dist[i][1])**2)**0.5)*(ela_dist[k+1][1]-ela_dist[k][1])/2*tw
        for j in range(i-1):
            if ela_dist[j+1][1]<=tf/d :
                area3=area3+((ela_dist[0][0]-(ela_dist[0][0]-ela_dist[0][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*ela_dist[i-1][1]*(ela_dist[j+1][1]-ela_dist[j][1])+0.5*(ela_dist[0][0]-ela_dist[0][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)*ela_dist[i-1][1]-0.5*(ela_dist[0][0]-(ela_dist[0][0]-ela_dist[0][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)-1/3*(ela_dist[0][0]-ela_dist[0][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**3-ela_dist[j][1]**3))*w
                #area3=area3+(ela_dist[j][0]*((ela_dist[j][1]-ela_dist[i-1][1])**2)**0.5+ela_dist[j+1][0]*((ela_dist[j+1][1]-ela_dist[i-1][1])**2)**0.5)*(ela_dist[j+1][1]-ela_dist[j][1])/2*w
            elif ela_dist[j][1]>=(d-tf)/d and ela_dist[j+1][1]<=1:
                area3=area3+((ela_dist[-1][0]-(ela_dist[-1][0]-ela_dist[-1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*ela_dist[i-1][1]*(ela_dist[j+1][1]-ela_dist[j][1])+0.5*(ela_dist[-1][0]-ela_dist[-1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)*ela_dist[i-1][1]-0.5*(ela_dist[-1][0]-(ela_dist[-1][0]-ela_dist[-1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)-1/3*(ela_dist[-1][0]-ela_dist[-1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**3-ela_dist[j][1]**3))*w
            elif ela_dist[j][1]>=tf/d and ela_dist[j+1][1]<=(d-tf)/d:
                area3=area3+((ela_dist[j][0]-(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*ela_dist[i-1][1]*(ela_dist[j+1][1]-ela_dist[j][1])+0.5*(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)*ela_dist[i-1][1]-0.5*(ela_dist[j][0]-(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)-1/3*(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**3-ela_dist[j][1]**3))*tw
                #area3=area3+(ela_dist[j][0]*((ela_dist[j][1]-ela_dist[i-1][1])**2)**0.5+ela_dist[j+1][0]*((ela_dist[j+1][1]-ela_dist[i-1][1])**2)**0.5)*(ela_dist[j+1][1]-ela_dist[j][1])/2*tw
        for k in range(i-1,len(ela_dist)-1):
            if ela_dist[k+1][1]<=tf/d :
                area4=area4+(0.5*(ela_dist[0][0]-(ela_dist[0][0]-ela_dist[0][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)+1/3*(ela_dist[0][0]-ela_dist[0][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**3-ela_dist[k][1]**3)-(ela_dist[0][0]-(ela_dist[0][0]-ela_dist[0][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*ela_dist[i-1][1]*(ela_dist[k+1][1]-ela_dist[k][1])-0.5*(ela_dist[0][0]-ela_dist[0][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)*ela_dist[i-1][1])*w 
                #area4=area4+(ela_dist[k][0]*((ela_dist[k][1]-ela_dist[i-1][1])**2)**0.5+ela_dist[k+1][0]*((ela_dist[k+1][1]-ela_dist[i-1][1])**2)**0.5)*(ela_dist[k+1][1]-ela_dist[k][1])/2*w
            elif ela_dist[k][1]>=(d-tf)/d and ela_dist[k+1][1]<=1:
                area4=area4+(0.5*(ela_dist[-1][0]-(ela_dist[-1][0]-ela_dist[-1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)+1/3*(ela_dist[-1][0]-ela_dist[-1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**3-ela_dist[k][1]**3)-(ela_dist[-1][0]-(ela_dist[-1][0]-ela_dist[-1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*ela_dist[i-1][1]*(ela_dist[k+1][1]-ela_dist[k][1])-0.5*(ela_dist[-1][0]-ela_dist[-1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)*ela_dist[i-1][1])*w 
            elif ela_dist[k][1]>=tf/d and ela_dist[k+1][1]<=(d-tf)/d:
                area4=area4+(0.5*(ela_dist[k][0]-(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)+1/3*(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**3-ela_dist[k][1]**3)-(ela_dist[k][0]-(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*ela_dist[i-1][1]*(ela_dist[k+1][1]-ela_dist[k][1])-0.5*(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)*ela_dist[i-1][1])*tw 
                #area4=area4+(ela_dist[k][0]*((ela_dist[k][1]-ela_dist[i-1][1])**2)**0.5+ela_dist[k+1][0]*((ela_dist[k+1][1]-ela_dist[i-1][1])**2)**0.5)*(ela_dist[k+1][1]-ela_dist[k][1])/2*tw    
        print('area1=',area1)
        print('area2=',area2)
        print('area3=',area3)
        print('area4=',area4)
        #area1=round(area1,11)
        #area2=round(area2,11)
        #area3=round(area3,11)
        #area4=round(area4,11)
        if area1>=area2 and area3<=area4:
            x = symbols('x',real=True)
            for j in range(i-1):
                if ela_dist[j+1][1]<=tf/d :
                    #area_t=area_t+(ela_dist[j][0]*((ela_dist[j][1]-x)**2)**0.5+ela_dist[j+1][0]*((ela_dist[j+1][1]-x)**2)**0.5)*(ela_dist[j+1][1]-ela_dist[j][1])/2*w
                    area_t=area_t+((ela_dist[0][0]-(ela_dist[0][0]-ela_dist[0][0])*(ela_dist[j][1]-ela_dist[j+1][1])**(-1)*ela_dist[j][1])*x*(ela_dist[j+1][1]-ela_dist[j][1])+0.5*(ela_dist[0][0]-ela_dist[0][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)*x-0.5*(ela_dist[0][0]-(ela_dist[0][0]-ela_dist[0][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)-1/3*(ela_dist[0][0]-ela_dist[0][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**3-ela_dist[j][1]**3))*w
                elif ela_dist[j][1]>=(d-tf)/d and ela_dist[j+1][1]<=1:
                    area_t=area_t+((ela_dist[-1][0]-(ela_dist[-1][0]-ela_dist[-1][0])*(ela_dist[j][1]-ela_dist[j+1][1])**(-1)*ela_dist[j][1])*x*(ela_dist[j+1][1]-ela_dist[j][1])+0.5*(ela_dist[-1][0]-ela_dist[-1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)*x-0.5*(ela_dist[-1][0]-(ela_dist[-1][0]-ela_dist[-1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)-1/3*(ela_dist[-1][0]-ela_dist[-1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**3-ela_dist[j][1]**3))*w
                elif ela_dist[j][1]>=tf/d and ela_dist[j+1][1]<=(d-tf)/d:
                    #area_t=area_t+(ela_dist[j][0]*((ela_dist[j][1]-x)**2)**0.5+ela_dist[j+1][0]*((ela_dist[j+1][1]-x)**2)**0.5)*(ela_dist[j+1][1]-ela_dist[j][1])/2*tw
                    area_t=area_t+((ela_dist[j][0]-(ela_dist[j][0]-ela_dist[j+1][0])*(ela_dist[j][1]-ela_dist[j+1][1])**(-1)*ela_dist[j][1])*x*(ela_dist[j+1][1]-ela_dist[j][1])+0.5*(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)*x-0.5*(ela_dist[j][0]-(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*ela_dist[j][1])*(ela_dist[j+1][1]**2-ela_dist[j][1]**2)-1/3*(ela_dist[j][0]-ela_dist[j+1][0])/(ela_dist[j][1]-ela_dist[j+1][1])*(ela_dist[j+1][1]**3-ela_dist[j][1]**3))*tw
            for k in range(i,len(ela_dist)-1):                                                                                                                       
                if ela_dist[k+1][1]<=tf/d :
                    area_b=area_b+(0.5*(ela_dist[0][0]-(ela_dist[0][0]-ela_dist[0][0])*(ela_dist[k][1]-ela_dist[k+1][1])**(-1)*ela_dist[k][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)+1/3*(ela_dist[0][0]-ela_dist[0][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**3-ela_dist[k][1]**3)-(ela_dist[0][0]-(ela_dist[0][0]-ela_dist[0][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*x*(ela_dist[k+1][1]-ela_dist[k][1])-0.5*(ela_dist[0][0]-ela_dist[0][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)*x)*w 
                    #area_b=area_b+(ela_dist[k][0]*((ela_dist[k][1]-x)**2)**0.5+ela_dist[k+1][0]*((ela_dist[j+1][1]-x)**2)**0.5)*(ela_dist[j+1][1]-ela_dist[j][1])/2*w
                elif ela_dist[k][1]>=(d-tf)/d and ela_dist[k+1][1]<=1:
                    area_b=area_b+(0.5*(ela_dist[-1][0]-(ela_dist[-1][0]-ela_dist[-1][0])*(ela_dist[k][1]-ela_dist[k+1][1])**(-1)*ela_dist[k][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)+1/3*(ela_dist[-1][0]-ela_dist[-1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**3-ela_dist[k][1]**3)-(ela_dist[-1][0]-(ela_dist[-1][0]-ela_dist[-1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*x*(ela_dist[k+1][1]-ela_dist[k][1])-0.5*(ela_dist[-1][0]-ela_dist[-1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)*x)*w 
                elif ela_dist[k][1]>=tf/d and ela_dist[k+1][1]<=(d-tf)/d:
                    area_b=area_b+(0.5*(ela_dist[k][0]-(ela_dist[k][0]-ela_dist[k+1][0])*(ela_dist[k][1]-ela_dist[k+1][1])**(-1)*ela_dist[k][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)+1/3*(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**3-ela_dist[k][1]**3)-(ela_dist[k][0]-(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*ela_dist[k][1])*x*(ela_dist[k+1][1]-ela_dist[k][1])-0.5*(ela_dist[k][0]-ela_dist[k+1][0])/(ela_dist[k][1]-ela_dist[k+1][1])*(ela_dist[k+1][1]**2-ela_dist[k][1]**2)*x)*tw 
                    #area_b=area_b+(ela_dist[j][0]*((ela_dist[j][1]-x)**2)**0.5+ela_dist[j+1][0]*((ela_dist[j+1][1]-x)**2)**0.5)*(ela_dist[j+1][1]-ela_dist[j][1])/2*tw
            print('area_t_o=',area_t)
            print('area_b_o=',area_b)
            print(ela_dist[i][1],ela_dist[i-1][1])
            if ela_dist[i][1] <= tf/d :
                area_t=area_t+((ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])*(ela_dist[i-1][1]-ela_dist[i][1])**(-1)*ela_dist[i-1][1])*x*(x-ela_dist[i-1][1])+0.5*(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*(x**2-ela_dist[i-1][1]**2)*x-0.5*(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*ela_dist[i-1][1])*(x**2-ela_dist[i-1][1]**2)-1/3*(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*(x**3-ela_dist[i-1][1]**3))*w
                #area_t=area_t+ela_dist[0][0]*((ela_dist[i-1][1]-x)**2)/2*w
                area_b=area_b+(0.5*(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])*(ela_dist[i-1][1]-ela_dist[i][1])**(-1)*ela_dist[i-1][1])*(ela_dist[i][1]**2-x**2)+1/3*(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*(ela_dist[i][1]**3-x**3)-(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*ela_dist[i-1][1])*x*(ela_dist[i][1]-x)-0.5*(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*(ela_dist[i][1]**2-x**2)*x)*w
            elif ela_dist[i-1][1]>=(d-tf)/d and ela_dist[i][1]<=1:
                area_t=area_t+((ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])*(ela_dist[i-1][1]-ela_dist[i][1])**(-1)*ela_dist[i-1][1])*x*(x-ela_dist[i-1][1])+0.5*(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*(x**2-ela_dist[i-1][1]**2)*x-0.5*(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*ela_dist[i-1][1])*(x**2-ela_dist[i-1][1]**2)-1/3*(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*(x**3-ela_dist[i-1][1]**3))*w
                #area_t=area_t+ela_dist[i-1][0]*((ela_dist[i-1][1]-x)**2)/2*w
                area_b=area_b+(0.5*(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])*(ela_dist[i-1][1]-ela_dist[i][1])**(-1)*ela_dist[i-1][1])*(ela_dist[i][1]**2-x**2)+1/3*(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*(ela_dist[i][1]**3-x**3)-(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*ela_dist[i-1][1])*x*(ela_dist[i][1]-x)-0.5*(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*(ela_dist[i][1]**2-x**2)*x)*w
            elif ela_dist[i-1][1]>=tf/d and ela_dist[i][1]<=(d-tf)/d:
                area_t=area_t+((ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])*(ela_dist[i-1][1]-ela_dist[i][1])**(-1)*ela_dist[i-1][1])*x*(x-ela_dist[i-1][1])+0.5*(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*(x**2-ela_dist[i-1][1]**2)*x-0.5*(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*ela_dist[i-1][1])*(x**2-ela_dist[i-1][1]**2)-1/3*(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*(x**3-ela_dist[i-1][1]**3))*tw
                area_b=area_b+(0.5*(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])*(ela_dist[i-1][1]-ela_dist[i][1])**(-1)*ela_dist[i-1][1])*(ela_dist[i][1]**2-x**2)+1/3*(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*(ela_dist[i][1]**3-x**3)-(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*ela_dist[i-1][1])*x*(ela_dist[i][1]-x)-0.5*(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]-ela_dist[i][1])*(ela_dist[i][1]**2-x**2)*x)*tw
            print('area_t=',area_t)
            print('area_b=',area_b)
            expr=area_t-area_b
            #print('Equation of neutral axis',expr)
            y_neu=solve(expr)
            y_neu_set=[]
            for sol in y_neu:
                sol=round(sol,14)
                y_neu_set.append(sol) 
            print('y_neu_set=',y_neu_set)  
            print('ela_dist2='+str(ela_dist[i][1])+', ela_dist1='+str(ela_dist[i-1][1]))
            for k in y_neu_set:
                if k<=ela_dist[i][1] and k>=ela_dist[i-1][1]: 
                    y_neutral=k
                    break
    print('y_neutral=',y_neutral)
    return ela_dist,y_neutral,tempdist
def E_Ix(ela_dist,y_neu,d,w,tf,tw,E):
    E_Ix=0
    E_Ix1=0
    for i in range(len(ela_dist)-1):
        if ela_dist[i+1][1] <= tf/d:
            E_Ix1=w*E*(0.25*(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*((ela_dist[i+1][1]*d)**4-(ela_dist[i][1]*d)**4)-1/3*(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*ela_dist[i][1]*d*((ela_dist[i+1][1]*d)**3-(ela_dist[i][1]*d)**3)+1/3*ela_dist[i][0]*((ela_dist[i+1][1]*d)**3-(ela_dist[i][1]*d)**3)-2/3*d*y_neu*(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*((ela_dist[i+1][1]*d)**3-(ela_dist[i][1]*d)**3)+(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*y_neu*d*ela_dist[i][1]*d*((ela_dist[i+1][1]*d)**2-(ela_dist[i][1]*d)**2)-ela_dist[i][0]*y_neu*d*((ela_dist[i+1][1]*d)**2-(ela_dist[i][1]*d)**2)+0.5*(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*((ela_dist[i+1][1]*d)**2-(ela_dist[i][1]*d)**2)*(y_neu*d)**2-(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*ela_dist[i][1]*d*(y_neu*d)**2*(ela_dist[i+1][1]*d-ela_dist[i][1]*d)+ela_dist[i][0]*(y_neu*d)**2*(ela_dist[i+1][1]-ela_dist[i][1])*d)
            E_Ix=E_Ix+E_Ix1
        elif ela_dist[i][1]>=(d-tf)/d:
            E_Ix1=w*E*(0.25*(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*((ela_dist[i+1][1]*d)**4-(ela_dist[i][1]*d)**4)-1/3*(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*ela_dist[i][1]*d*((ela_dist[i+1][1]*d)**3-(ela_dist[i][1]*d)**3)+1/3*ela_dist[i][0]*((ela_dist[i+1][1]*d)**3-(ela_dist[i][1]*d)**3)-2/3*d*y_neu*(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*((ela_dist[i+1][1]*d)**3-(ela_dist[i][1]*d)**3)+(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*y_neu*d*ela_dist[i][1]*d*((ela_dist[i+1][1]*d)**2-(ela_dist[i][1]*d)**2)-ela_dist[i][0]*y_neu*d*((ela_dist[i+1][1]*d)**2-(ela_dist[i][1]*d)**2)+0.5*(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*((ela_dist[i+1][1]*d)**2-(ela_dist[i][1]*d)**2)*(y_neu*d)**2-(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*ela_dist[i][1]*d*(y_neu*d)**2*(ela_dist[i+1][1]*d-ela_dist[i][1]*d)+ela_dist[i][0]*(y_neu*d)**2*(ela_dist[i+1][1]-ela_dist[i][1])*d)
            E_Ix=E_Ix+E_Ix1   
        elif ela_dist[i][1] >= tf/d and ela_dist[i+1][1]<=(d-tf)/d:
            E_Ix1=tw*E*(0.25*(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*((ela_dist[i+1][1]*d)**4-(ela_dist[i][1]*d)**4)-1/3*(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*ela_dist[i][1]*d*((ela_dist[i+1][1]*d)**3-(ela_dist[i][1]*d)**3)+1/3*ela_dist[i][0]*((ela_dist[i+1][1]*d)**3-(ela_dist[i][1]*d)**3)-2/3*d*y_neu*(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*((ela_dist[i+1][1]*d)**3-(ela_dist[i][1]*d)**3)+(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*y_neu*d*ela_dist[i][1]*d*((ela_dist[i+1][1]*d)**2-(ela_dist[i][1]*d)**2)-ela_dist[i][0]*y_neu*d*((ela_dist[i+1][1]*d)**2-(ela_dist[i][1]*d)**2)+0.5*(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*((ela_dist[i+1][1]*d)**2-(ela_dist[i][1]*d)**2)*(y_neu*d)**2-(ela_dist[i][0]-ela_dist[i+1][0])/((ela_dist[i][1]*d-ela_dist[i+1][1]*d))*ela_dist[i][1]*d*(y_neu*d)**2*(ela_dist[i+1][1]*d-ela_dist[i][1]*d)+ela_dist[i][0]*(y_neu*d)**2*(ela_dist[i+1][1]-ela_dist[i][1])*d)
            E_Ix=E_Ix+E_Ix1
    return E_Ix
def cen_of_mass(a,b,h):
    a1=max(a,b)
    b1=min(a,b)
    c_mass=(a1+2*b1)*h/3/(a1+b1)
    return c_mass
def deformation_nonuniform(temp_dist,ela_dist,alpha,d,w,tf,tw,ela_axis):
    eps_dist=[]
    for i in range(len(temp_dist)):
        eps_dist.append([(temp_dist[i][0]-20)*alpha,temp_dist[i][1]])
    for i in range(len(eps_dist)-1):
        if eps_dist[i][1]<=ela_axis and eps_dist[i+1][1]>ela_axis:
            eps_ela=interpolation(ela_axis,eps_dist)
            eps_dist.insert(i+1,[eps_ela,ela_axis])
            break
    print('eps_dist=',eps_dist)
    F=0
    M=0
    b = symbols('b',real=True)
    k = symbols('k',real=True)
    print('len(ela_dist)=',len(ela_dist))
    print('len(eps_dist)=',len(eps_dist))
    for i in range(1,len(ela_dist)):
        if ela_dist[i][1]<=tf/d or ela_dist[i-1][1]>=1-tf/d:
            F=F+(1/3*(ela_dist[i][0]-ela_dist[i-1][0])/(ela_dist[i][1]*d-ela_dist[i-1][1]*d)*((eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)-k)*((ela_dist[i][1]*d)**3-(ela_dist[i-1][1]*d)**3)+0.5*(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]*d-ela_dist[i][1]*d)*ela_dist[i-1][1]*d)*((eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)-k)*((eps_dist[i][1]*d)**2-(eps_dist[i-1][1]*d)**2)+0.5*(ela_dist[i][0]-ela_dist[i-1][0])/(ela_dist[i][1]*d-ela_dist[i-1][1]*d)*(eps_dist[i-1][0]-(eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)*eps_dist[i-1][1]*d-b)*((eps_dist[i][1]*d)**2-(eps_dist[i-1][1]*d)**2)+(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]*d-ela_dist[i][1]*d)*ela_dist[i-1][1]*d)*(eps_dist[i-1][0]-(eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)*ela_dist[i-1][1]*d-b)*(ela_dist[i][1]*d-ela_dist[i-1][1]*d))*w
            M=M+(1/4*(ela_dist[i][0]-ela_dist[i-1][0])/(ela_dist[i][1]*d-ela_dist[i-1][1]*d)*((eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)-k)*((ela_dist[i][1]*d)**4-(ela_dist[i-1][1]*d)**4)+1/3*(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]*d-ela_dist[i][1]*d)*ela_dist[i-1][1]*d)*((eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)-k)*((eps_dist[i][1]*d)**3-(eps_dist[i-1][1]*d)**3)+1/3*(ela_dist[i][0]-ela_dist[i-1][0])/(ela_dist[i][1]*d-ela_dist[i-1][1]*d)*(eps_dist[i-1][0]-(eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)*eps_dist[i-1][1]*d-b)*((eps_dist[i][1]*d)**3-(eps_dist[i-1][1]*d)**3)+0.5*(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]*d-ela_dist[i][1]*d)*ela_dist[i-1][1]*d)*(eps_dist[i-1][0]-(eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)*eps_dist[i-1][1]*d-b)*((ela_dist[i][1]*d)**2-(ela_dist[i-1][1]*d)**2))*w        
        elif ela_dist[i-1][1]>tf/d and ela_dist[i][1]<=(d-tf)/d:
            F=F+(1/3*(ela_dist[i][0]-ela_dist[i-1][0])/(ela_dist[i][1]*d-ela_dist[i-1][1]*d)*((eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)-k)*((ela_dist[i][1]*d)**3-(ela_dist[i-1][1]*d)**3)+0.5*(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]*d-ela_dist[i][1]*d)*ela_dist[i-1][1]*d)*((eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)-k)*((eps_dist[i][1]*d)**2-(eps_dist[i-1][1]*d)**2)+0.5*(ela_dist[i][0]-ela_dist[i-1][0])/(ela_dist[i][1]*d-ela_dist[i-1][1]*d)*(eps_dist[i-1][0]-(eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)*eps_dist[i-1][1]*d-b)*((eps_dist[i][1]*d)**2-(eps_dist[i-1][1]*d)**2)+(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]*d-ela_dist[i][1]*d)*ela_dist[i-1][1]*d)*(eps_dist[i-1][0]-(eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)*ela_dist[i-1][1]*d-b)*(ela_dist[i][1]*d-ela_dist[i-1][1]*d))*tw
            M=M+(1/4*(ela_dist[i][0]-ela_dist[i-1][0])/(ela_dist[i][1]*d-ela_dist[i-1][1]*d)*((eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)-k)*((ela_dist[i][1]*d)**4-(ela_dist[i-1][1]*d)**4)+1/3*(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]*d-ela_dist[i][1]*d)*ela_dist[i-1][1]*d)*((eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)-k)*((eps_dist[i][1]*d)**3-(eps_dist[i-1][1]*d)**3)+1/3*(ela_dist[i][0]-ela_dist[i-1][0])/(ela_dist[i][1]*d-ela_dist[i-1][1]*d)*(eps_dist[i-1][0]-(eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)*eps_dist[i-1][1]*d-b)*((eps_dist[i][1]*d)**3-(eps_dist[i-1][1]*d)**3)+0.5*(ela_dist[i-1][0]-(ela_dist[i-1][0]-ela_dist[i][0])/(ela_dist[i-1][1]*d-ela_dist[i][1]*d)*ela_dist[i-1][1]*d)*(eps_dist[i-1][0]-(eps_dist[i-1][0]-eps_dist[i][0])/(eps_dist[i-1][1]*d-eps_dist[i][1]*d)*eps_dist[i-1][1]*d-b)*((ela_dist[i][1]*d)**2-(ela_dist[i-1][1]*d)**2))*tw   
    var=solve([F,M])
    b=var[b]
    k=var[k]
    return b,k,eps_dist
def thermal_bowing(b,k,d,L,x,BC):    
    rou=abs((1+b+k*0.5*d)*d/(k*d))
    print('rou=',rou)
    if BC==0:
        theta=abs(L*(k*d)/d)
    elif BC==1:
        theta=0
        rou=9999999999999
    u1=((rou)**2-(rou*math.sin(theta/2)-x)**2)**0.5
    u=u1-(rou*math.cos(theta/2))
    return u
def range_x(b,k,d,L):
    rou=abs((1+(b+k*0.5*d))*d/(k*d))
    theta=abs(L*(k*d)/d)
    X=2*(rou)*math.sin(theta/2)
    return X
    
def M_bowing(ela_dist,ela_axis,d,w,tf,tw,E):
    M_b=0
    for i in range(len(ela_dist)-1):
        if ela_dist[i+1][1]<=tf/d or ela_dist[i][1]>=1-tf/d:
            h=ela_dist[i+1][1]-ela_dist[i][1]
            arm=ela_axis-(h-cen_of_mass(ela_dist[i][0],ela_dist[i+1][0],h))
            M_b=M_b+0.5*(ela_dist[i][0]+ela_dist[i+1][0])*h*arm*w*d**2*E
        elif ela_dist[i][1]>tf/d and ela_dist[i+1][1]<=(d-tf)/d:
            h=ela_dist[i+1][1]-ela_dist[i][1]
            arm=ela_axis-(h-cen_of_mass(ela_dist[i][0],ela_dist[i+1][0],h))
            M_b=M_b+0.5*(ela_dist[i][0]+ela_dist[i+1][0])*h*arm*tw*d**2*E    
    return M_b
def cal_u(M_b,x,L,E_Ix):
    u_pp=M_b/E_Ix
    print('u_pp=',u_pp)
    u=0.5*u_pp*x**2-L/2*u_pp*x
    return u
    
n=0
chi_LT_fi_11=[]
chi_LT_fi_8=[]
chi_LT_fi_6=[]
chi_LT_fi_35=[]
M_cr_group=[]
M_cr_cluster=[]
lamda4=[]
lamda3=[]
lamda2=[]
lamda1=[]
kk=0
ela_axis_set=[]
#def M_b(length,temp,a,b,c,h,n):
#for L in length:
mpl_group=[]
lamda_group=[]
chi_LT_fi_group=[]
d=258.3
w=146.1
tf=9.1
tw=6.1
a1=tf/d
a2=tf/d+0.00001
a3=1-tf/d-0.00001
a4=1-tf/d
for tp in temp:
    M_cr_set=[]
    tempdist_set=[]
    tempdist_set.append([[tp,0],[tp,0.5],[tp,1]])
    tempdist_set.append([[tp,0],[tp*0.75,0.5],[tp*0.5,1]])
    #tempdist_set.append([[tp,0],[tp,a1],[tp*0.75,a2],[tp*0.75,a3],[tp*0.5,a4],[tp*0.5,1]])
    tempdist_set.append([[tp,0],[tp*0.4,0.5],[tp*0.2,1]])
    tempdist_set.append([[tp,0],[tp*0.2,0.5],[tp*0.1,1]])
    lamda_set=[]
    M_b_fi_set=[]
    mpl_set=[]
    chi_LT_fi_set=[]
    for j in range(1,2):
        tempdist=tempdist_set[j]
        tempdist=[[1141.39,0],[940.532,0.194735],[753.094,0.353465],[614.9,0.516067],[511.181,0.678668],[429.928,0.837398],[365.723,1]]
        #tempdist=[[500,0],[200,0.5],[100,1]]
        tempdist=[[500,0],[100,0.5],[50,1]]
        #stre_upper=stre_red(tempdist[0][0]*tempdist[1][0])
        Iz=tf*(w**3)/12*2+(d-2*tf)*(tw**3)/12
        Iw=Iz*(d-tf)**2/4   
        #It=1.578e5
        It=(2*w*tf**3+(d-tf)*tw**3)/3
        print('It='+str(It))
        #Wy=1e5
        fpl=320
        E=210000
        G=80770
        C1=1.6
        C2=0.454
        zg=0
        alpha=1e-5
        k=1
        kw=1
        Wpl=w*tf*(d-tf)+0.25*tw*(d-2*tf)**2
        # zg=
        gamma_M1=''
        tempdist_uf=[]
        tempdist_lf=[]
        tempdist_web=[]
        tempdist_uf.append(tempdist[0])
        tempdist_uf.append([tempdist[0][0],1])
        tempdist_lf.append([tempdist[-1][0],0])
        tempdist_lf.append(tempdist[-1])
        tempdist_web.append([tempdist[1][0],0])
        tempdist_web.append([tempdist[1][0],1])
        #print(tempdist)
        ela_dist,ela_axis,tempdist=elastic_neutral_axis(d,w,tf,tw,tempdist,fpl,E,Iz,Iw,It,Wpl,G)
        ela_axis_set.append(ela_axis)
        for i in range(len(ela_dist)-1):
            if ela_dist[i][1]<=ela_axis and ela_dist[i+1][1]>ela_axis:
                k_ela=interpolation(ela_axis,ela_dist)
                ela_dist.insert(i+1,[k_ela,ela_axis])
                break
        b,k,eps_dist=deformation_nonuniform(tempdist,ela_dist,alpha,d,w,tf,tw,ela_axis)
        deform_eps=[]
        print('b='+str(b)+', k='+str(k))
        y_range=np.linspace(0,1,100)
        for y in y_range:
            deform_eps.append(k*d*y+b)
        deform_eps=np.array(deform_eps)
        plt.plot(deform_eps,y_range,label='Deformed section')
        eps_range=[]
        y_eps_range=[]
        for x in eps_dist:
            eps_range.append(x[0])
            y_eps_range.append(x[1])
        plt.plot(eps_range,y_eps_range,label='Undeformed section',linestyle='--')
        plt.ylim(1,0)
        plt.ylabel('Deepness from top/depth y/d',fontsize=12,fontstyle='italic',fontweight='bold')
        plt.xlabel('Strain ε',fontsize=12,fontstyle='italic',fontweight='bold')
        plt.legend()
        plt.savefig('deformed plane.png')
        plt.show()
        print('b='+str(b)+', k='+str(k))
        X_range=range_x(b,k,d,6000)
        X_range=float(X_range)
        print('X_range=',X_range)
        x_array=np.linspace(0,X_range,100)
        u_array=[]
        for x in x_array:
            u=thermal_bowing(b,k,d,6000,x,0)
            u_array.append(u)
            u_max=max(u_array)
        print('u_max=',u_max)
        u_array=np.array(u_array)
        with open('thermal_bowing'+str(j)+'.csv','a') as f:
            for i in range(len(x_array)):
                f.write(str(x_array[i])+','+str(u_array[i])+'\n')
        plt.plot(x_array,u_array)
        plt.show()
        
print(ela_axis_set)