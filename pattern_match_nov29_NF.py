# -*- coding: utf-8 -*-
#Will Cole, 2016, will.cole3@gmail.com

#Updated last: Nov 29, 2016 
#Added function to match adjacent gsK sets
#function is called matchK
#altered main() to included both new and old versions
#changed the tolerances for P/Q/R matches

#This is a pattern matching program designed for the D2O dimer, though is capable of working for any system so long as you change the constants
#It currently does not incorporate tunneling stucture, that is done by hand (2016)
#It takes an input of a line list in form [Wavenumber, intensity]
#It will split out the possible matches for P/Q/R branches that it finds


###################Instructions#################################################
#You can run this two ways: in a distribution or from terminal
    #If you run in dist. the program will ask for constants in the dist.
    #terminal.
        
    #If you run from terminal the program will ask the user for the values of the constants
        #the call will look like: python pattern_match.py 
        #see below for constant definitions; maxgsKvalue is the max value of K in the ground state (normally 2-10)  
import time
import numpy as np
import datetime as ds

today=ds.date.today()

#This is the raw line list
#The format must be [line (W#),S/N]
#data=np.loadtxt('remaininglines.txt')

#Don't touch these variables here
global upper2B 
global lower2B
global wiggle
global lowIntTol
global hiIntTol
global diffTol
global inM
global label
#THIS IS FOR TESTING
#DONT TOUCH
################################################################################
#NO TOUCH(all values should be in the same units as those the lines in data [W# or Hz])
#Now is a location for some constants

#These are the upper and lower values for 2B
#this determines the system you want to look at
lower2B = 0
upper2B = 0

#This is the tolerance for matching, typical value is 30-50MHz
wiggle = 0

#These are the tolerances for the intensity ratios of adjacent lines
#this is based on the getLineInts func below and basically represents how 
#confident in those predicted intensities
lowIntTol=0
hiIntTol=0

#This is a tolerance for the delB values in the Q branches
diffTol=0

################################################################################
#STOP TOUCHING, go to line 800ish for main func def

def indexData(arr):
    """Give each data row a number, this is used for keeping track of each line
        
    Parameters
    ----------
    arr : the array to index
    
    Returns
    -------
    data : the input array, arr, with a new column consisting of row number
    
    """
    size=len(arr)
    data=np.zeros([size,3])
    for i in range(len(arr)):
        data[i,0]=arr[i,0]
        data[i,1]=arr[i,1]
        data[i,2]=i
    return data



#Want something that can take the first line in a P/Q/R branch gen a intensity pattern
#There is some choice here, in where you think the intensity maximum will occur in a PQR stack
def getLineInts(j, size):
    """Get a line intensity pattern accounting for degeneracy and a Boltzmann temp profile:
    
    Parameters
    ---------
    j : the j value of the transition you expect to be highest in intensity; this is essentially
        a guess, you can hedge by using known ground state transition energies and assuming a
        rotational temp of ~ 4K
    
    size : the largest j value to calculate
    
    Returns
    -------
    intArray : a 2D array with format: j-value, relative intensity
        
    """
    intArray=np.zeros([size,2])
    first=((2*j)+1)*np.exp((-1*j)/6.1)
    for i in range(0,size):
        intArray[i,0]=i
        #use normal deg factor with exp decay after J=critj
        #forced to be normalized to the first line
        lineInt=(((2*(j+i))+1)*np.exp((-1*(j+i))/6.1))/first
        intArray[i,1]=lineInt
    return intArray

#Literally a kill function to remove line after matched if we want to do that
#That ends up being a very bad idea
def killLine(orgArr,matchArr):
    """DO NOT use this function right now, might use later Apr. 2016
        
    Parameters
    ----------
    orgArr : the original array to remove lines from
    
    matchArr : the matches found by the program, these are the lines to remove from orgArr
    
    Returns
    -------
    remain : the original array with the math transitions removed
        
    """
    remain=np.zeros([0,3])
    for i in range(len(orgArr)):
        #for j in range(len(matchArr)):
            if orgArr[i,2] not in matchArr[:,3]:
                
            
                dummy=np.zeros([1,3])
                
                for k in range(0,3):
                    dummy[0,k]=orgArr[i,k]
                remain=np.vstack([remain,dummy])
    return remain


#default func for writing out a array to file
def writeOut(self, inArr, end):
    """Write out the list of possible match found by the algorithm
        
    Parameters
    ---------
    self : name of the output file
    
    inArr : the array containing the possible matches, formating follows the header below
    
    end : width of the inArr
        
    """
    f=open(self+'_'+'%s' %today+'.txt','w')
    f.write('The column headings are as follows:'+'\n')
    f.write('"Line (W#)" is the line in units of wavenumber'+'\n')
    f.write('S/N is line signal to noise'+'\n')
    f.write('BS distinguishes different possible P/Q/R sets'+'\n')
    f.write('ODR is "Original Data Row" helps you find the line in raw data' +'\n')
    f.write('J Val is the g.s. J value of the line'+'\n')
    f.write('Diff is 2B for P/R or 2delB for Q branches' +'\n')
    f.write('P/Q/R tells you what type each branch set is -1=P, 0=Q, 1=R' +'\n')
    f.write('gsK is the ground state K value'+'\n')
    f.write('Set label is a label for the stitchtog output'+'\n')
    f.write('Tol is the delM value found'+'\n')
    #This is the default format for the headers
    format=['Line (W#)', 'S/N', 'BS', 'ODR', 'J Val','Diff', 'P/Q/R','gsK','Set Label','Tol.']
    f.write('Output file for '+self+'\n')
    #First want to give headers to all columns
    for i in range(len(format)):
        f.write(format[i]+'\t')
    f.write('\n')
    #Now we are ready to begin writing the actual data from the inArr
    for j in range(len(inArr)):
        
        if inArr[j,0] == -1:
            f.write('New Set')
            f.write('\n')
        else:
            for k in range(0,end):
                f.write('%s' %inArr[j,k]+'\t')
            f.write('\n')
    f.close()

#want a function that searches a given array to find possible P/Q/R branches
#the inArr needs to have the freq on in 0th col and intensity on 1st col
#The outArr will be a m x 8 array 
#Have a label for the gsK that will be fed in by Stitchtog
def findRBranch(inArr,startj, tol,gsK):
    """Method to find transitions which met the criteria for the an R branch
        
    Parameters
    ----------
    inArr : the list of transitions to test
    
    startj : the starting j value
    
    tol : the tolerance for the match, user defined
    
    gsK : the ground state K value
    
    Returns
    -------
    outArr : the possible matchs formated in the following way: transition, transition intensity,
        match number, orgArr row number, start j value, R indicator (=1), ground state K value
    
    inArr : the original input array
        
    """
    outArr=np.zeros([0,8])
    tab=0
    refIntArray=getLineInts(startj,10)
    for i in range(len(inArr)):
        dummy=np.zeros([1,8])
        first=inArr[i,0]
        refInt=inArr[i,1]
        dummy[0,0]=inArr[i,0]
        dummy[0,1]=inArr[i,1]
        dummy[0,2]=tab
        dummy[0,3]=inArr[i,2]
        dummy[0,4]=startj
        dummy[0,6]=1
        dummy[0,7]=gsK
        for j in range(len(inArr)):
            #These ranges are set by the known 2B spacing of water clusters (dimer to octamer)
            if i < j:
                second=inArr[j,0]
                testint=(inArr[j,1]/refInt)
                if lower2B <= (second-first) <= upper2B and (lowIntTol*refIntArray[1,1])<=testint<=(hiIntTol*refIntArray[1,1]):
                    dummy2=np.zeros([1,8])
                    dummy2[0,0]=second
                    dummy2[0,1]=inArr[j,1]
                    dummy2[0,2]=tab
                    dummy2[0,3]=inArr[j,2]
                    dummy2[0,4]=startj+1
                    dummy2[0,6]=1
                    dummy2[0,7]=gsK
                    diff=second-first
                    dummy[0,5]=round(diff,4)
                    dummy2[0,5]=round(diff,4)
                    run=np.zeros([0,8])
                    for k in range(1,8):
                        for m in range(len(inArr)):
                            if j < m:
                                third=inArr[m,0]
                                secdiff=((third-second))
                                testint=(inArr[m,1]/refInt)
                                addOn=(((tol)*(k**2))/diff)
                                if (k-addOn) <= (secdiff/(diff+0.00000001)) <= (k+addOn) and (lowIntTol*refIntArray[(1+k),1])<=testint<=(hiIntTol*refIntArray[(1+k),1]):
                                
                                    dummy3=np.zeros([1,8])
                                    dummy3[0,0]=third
                                    dummy3[0,1]=inArr[m,1]
                                    dummy3[0,2]=tab
                                    dummy3[0,3]=inArr[m,2]
                                    dummy3[0,4]=(startj+1)+k
                                    dummy3[0,5]=round((secdiff/k),4)
                                    dummy3[0,6]=1
                                    dummy3[0,7]=gsK
                                    run=np.vstack([run,dummy3])
                                    j=m
                                    break
                    if 1 <= len(run):
                        outArr=np.vstack([outArr,dummy])
                        outArr=np.vstack([outArr,dummy2])
                        outArr=np.vstack([outArr,run])
                        tab+=1
                        break      
    return outArr, inArr

#These are test points for each branch                               

#out,inn=findRBranch(data,0,wiggle)  
#print len(out)         
                    
                
def findPBranch(inArr,startj, tol,gsK):
    """Method to find transitions which met the criteria for the an P branch
        
        Parameters
        ----------
        inArr : the list of transitions to test
        
        startj : the starting j value
        
        tol : the tolerance for the match, user defined
        
        gsK : the ground state K value
        
        Returns
        -------
        outArr : the possible matchs formated in the following way: transition, transition intensity,
        match number, orgArr row number, start j value, P indicator (= -1), ground state K value
        
        inArr : the original input array
        
        """
    outArr=np.zeros([0,8])
    tab=0
    refIntArray=getLineInts(startj,10)
    for i in range(len(inArr)):
        dummy=np.zeros([1,8])
        first=inArr[i,0]
        refInt=inArr[i,1]
        dummy[0,0]=inArr[i,0]
        dummy[0,1]=inArr[i,1]
        dummy[0,2]=tab
        dummy[0,3]=inArr[i,2]
        dummy[0,4]=startj
        dummy[0,6]=-1
        dummy[0,7]=gsK
        for j in range(len(inArr)):
            #These ranges are set by the known 2B spacing of water clusters (dimer to octamer)
            if j < i:
                second=inArr[j,0]
                testint=(inArr[j,1]/refInt)
                if lower2B <= (first-second) <= upper2B and (lowIntTol*refIntArray[1,1])<=testint<=(hiIntTol*refIntArray[1,1]):
                
                    dummy2=np.zeros([1,8])
                    dummy2[0,0]=second
                    dummy2[0,1]=inArr[j,1]
                    dummy2[0,2]=tab
                    dummy2[0,3]=inArr[j,2]
                    dummy2[0,4]=(1+startj)
                    dummy2[0,6]=-1
                    dummy2[0,7]=gsK
                    diff=first-second
                    dummy[0,5]=round(diff,4)
                    dummy2[0,5]=round(diff,4)
                    run=np.zeros([0,8])
                    for k in range(1,8):
                        for m in range(len(inArr)):
                            if m < j:
                                third=inArr[m,0]
                                secdiff=((second-third))
                                testint=(inArr[m,1]/refInt)
                                addOn=(((tol)*(k**2))/diff)
                                if (k-addOn) <= (secdiff/(diff+0.00000001)) <= (k+addOn) and (lowIntTol*refIntArray[(1+k),1])<=testint<=(hiIntTol*refIntArray[(1+k),1]):
                                
                                    dummy3=np.zeros([1,8])
                                    dummy3[0,0]=third
                                    dummy3[0,1]=inArr[m,1]
                                    dummy3[0,2]=tab
                                    dummy3[0,3]=inArr[m,2]
                                    dummy3[0,4]=(1+startj)+k
                                    dummy3[0,5]=round((secdiff/k),4)
                                    dummy3[0,6]=-1
                                    dummy3[0,7]=gsK
                                    run=np.vstack([run,dummy3])
                                    j=m
                                    break
                    if 1 <= len(run):
                        outArr=np.vstack([outArr,dummy])
                        outArr=np.vstack([outArr,dummy2])
                        outArr=np.vstack([outArr,run])
                        tab+=1
                        break            
    return outArr, inArr                            
            
#Test point for P branch     

#out,inn=findPBranch(data,1,wiggle)  
#print len(out)

        
def findQBranch(inArr,startj, tol,gsK):
    """Method to find transitions which met the criteria for the an Q branch
        
        Parameters
        ----------
        inArr : the list of transitions to test
        
        startj : the starting j value
        
        tol : the tolerance for the match, user defined
        
        gsK : the ground state K value
        
        Returns
        -------
        outArr : the possible matchs formated in the following way: transition, transition intensity,
        match number, orgArr row number, start j value, Q indicator (=0), ground state K value
        
        inArr : the original input array
        
        """
    outArr=np.zeros([0,8])
    tab=0
    refIntArray=getLineInts(startj,10)
    for i in range(len(inArr)):
        dummy=np.zeros([1,8])
        first=inArr[i,0]
        refInt=inArr[i,1]
        dummy[0,0]=inArr[i,0]
        dummy[0,1]=inArr[i,1]
        dummy[0,2]=tab
        dummy[0,3]=inArr[i,2]
        dummy[0,4]=startj
        dummy[0,6]=0
        dummy[0,7]=gsK
        for j in range(len(inArr)):
            #These ranges are set by the known 2B spacing of water clusters (dimer to octamer)
            if j != i:
                second=inArr[j,0]
                testint=(inArr[j,1]/refInt)
                if (0.95*((2*tol*startj)+(2*tol))) <= (second-first) <= (1.05*((2*tol*startj)+(2*tol))) or -1*(1.05*((2*tol*startj)+(2*tol))) <= (second-first) <= -1*(0.95*((2*tol*startj)+(2*tol))) and ((lowIntTol)*refIntArray[1,1])<=testint<=((hiIntTol)*refIntArray[1,1]):
                
                    dummy2=np.zeros([1,8])
                    dummy2[0,0]=second
                    dummy2[0,1]=inArr[j,1]
                    dummy2[0,2]=tab
                    dummy2[0,3]=inArr[j,2]
                    dummy2[0,4]=startj+1
                    dummy2[0,6]=0
                    dummy2[0,7]=gsK
                    diff=second-first
                    dummy[0,5]=round(diff,4)
                    dummy2[0,5]=round(diff,4)
                    run=np.zeros([0,8])
                    for k in range(1,8):
                        for m in range(len(inArr)):
                            if m != j:
                                secJ=dummy2[0,4]
                                divfact=((k*(secJ+((k+1)/2)))/(secJ+1))
                                third=inArr[m,0]
                                secdiff=((third-second)/divfact)
                                testint=(inArr[m,1]/refInt)
                                #the 0.00000001 is just to avoid dividing by zero
                                if 0.9 <= (secdiff/(diff+0.00000001)) <= 1.1 and ((lowIntTol)*refIntArray[(1+k),1])<=testint<=((hiIntTol)*refIntArray[(1+k),1]):
                                
                                    dummy3=np.zeros([1,8])
                                    dummy3[0,0]=third
                                    dummy3[0,1]=inArr[m,1]
                                    dummy3[0,2]=tab
                                    dummy3[0,3]=inArr[m,2]
                                    dummy3[0,4]=(startj+1)+k
                                    dummy3[0,5]=round(secdiff,4)
                                    dummy3[0,6]=0
                                    dummy3[0,7]=gsK
                                    run=np.vstack([run,dummy3])
                                    j=m
                                    break
                    if 1 <= len(run):
                        outArr=np.vstack([outArr,dummy])
                        outArr=np.vstack([outArr,dummy2])
                        outArr=np.vstack([outArr,run])
                        tab+=1
                        break       
    return outArr, inArr                            
            
#test point for Q branch       

#out,inn=findQBranch(data,0,wiggle)  
#writeOut('test',out)


def prBSpacing(B, rJ,pJ):
    """Get the expected between r and p branched for given j values
        
    Parameters
    ----------
    B : The given B value
    
    rJ : the J value of the first r branch transition
    
    pJ : the J value for the first p branch transition
    
    Returns
    -------
    outLen : the expected spacing between P and R branches
        
    """
    ref2B=B
    space=1+(rJ+pJ)
    outLen=(space)*ref2B
    return outLen

def qBSpacing(B, rJ):
    """Get the expected between r and q branched for given j values
        
        Parameters
        ----------
        B : The given B value
        
        rJ : the J value of the first r branch transition
        
        
        Returns
        -------
        outLen : the expected spacing between Q and R branches
        
        """
    ref2B=B
    space=1+(rJ)
    outLen=(space)*ref2B
    return outLen   

#This func takes a matched set and writes the branch set to an outArr
#Each possible match will be seperated by line of -1's
def zipSet(inRArr, firstRRow, inPArr, firstPRow, inQArr, firstQRow):
    """Method to mesh matched R, P, and Q branches
        
    Parameters
    ----------
    inRArr : The matched R branch transitions
    
    firstRRow : the index of the first R transition
    
    inPArr : The matched P branch transitions
    
    firstProw : the index of the first P transition
    
    inQArr : The matched Q branch transitions
    
    firstQRow : the index of the first Q transition
    
    Returns
    -------
    outArr: The matched R,P, and Q branches
        
    """
    outArr=np.zeros([0,8])
    if len(inQArr) == 0 and firstQRow == -1:
        breaker=np.array([-1,-1,-1,-1,-1,-1,-1,-1])
        outArr=np.vstack([outArr,breaker])
        for i in range(firstRRow, len(inRArr)):
            
            if inRArr[i,2] == inRArr[firstRRow,2]:
                dummy=inRArr[i,:]
                outArr=np.vstack([outArr, dummy])
        for i in range(firstPRow, len(inPArr)):
            
            if inPArr[i,2] == inPArr[firstPRow,2]:
                dummy=inPArr[i,:]
                outArr=np.vstack([outArr, dummy])  
    else:
        breaker=np.array([-1,-1,-1,-1,-1,-1,-1,-1])
        outArr=np.vstack([outArr,breaker])
        for i in range(firstRRow, len(inRArr)):
            
            if inRArr[i,2] == inRArr[firstRRow,2]:
                dummy=inRArr[i,:]
                outArr=np.vstack([outArr, dummy])
        for i in range(firstQRow, len(inQArr)):
            
            if inQArr[i,2] == inQArr[firstQRow,2]:
                dummy=inQArr[i,:]
                outArr=np.vstack([outArr, dummy])
        for i in range(firstPRow, len(inPArr)):
            
            if inPArr[i,2] == inPArr[firstPRow,2]:
                dummy=inPArr[i,:]
                outArr=np.vstack([outArr, dummy])
    return outArr


#Now that we have a way to search for all the possible P/Q/R branches
#We need another func to gen and search through those array and piece
#together each P with a Q and R subject to the selection rules imposed

#The normal selection rules are:
#Parallel type-> delK=0 in two flavors
    #K"=0 -> no Q branch
    #K"=/=0 -> has Q branch
#perpendicular type -> delK=±1, both have Q branches

#each of these types has a unique relationship between the first line in the 
#P/Q/R branches, which can be used to assign the type AND K" which is all we need to 
#make a definitive assignment.

#Also is would probabaly be benificial to work system by system to avoid removing potential transitions from
#the most likely systems.

#This will gen the appropriate list of possible transitions and then look for 
#full P/Q/R branches that obey selection rules
#Need to update this to have a label for the gsK in the outArr

def stitchTog(inArr,delK, gsK, tol):
    """Function to text R, P, and Q branches agaisnt each other to determine if they belong to a common set
        
    Parameters
    ---------
    inArr : the array of transitions found experimentally
    
    delK : the delta K values to test (can be -1, 0, or 1)
    
    gsK : the ground state K value
    
    tol : the user defined tolerance for the match
    
    Returns
    -------
    outArr : the array of matched R, P, and Q transitions
    
    inArr : the original array inArr
        
    """
    global label
    outArr=np.zeros([0,9])
    #The first step is to find the appropriate case
    #make sure you don't violate selection rules
    if delK == 0:
        #special case for delK=0, K"=0
        if gsK == 0:
            posr,inr=findRBranch(inArr, gsK, wiggle,gsK)
            posp,inp=findPBranch(inArr, (gsK+1),wiggle,gsK)
            for i in range(len(posr)):
                if posr[i,4] == gsK:
                    firstr=posr[i,0]
                    firstrint=posr[i,1]
                    #this was the matched 2B value
                    testDiff=posr[i,5]
                    for j in range(len(posp)):
                        if posp[j,4] == (gsK+1) and 0.975 <= (testDiff/posp[j,5]) <= 1.025:
                            firstp=posp[j,0]
                            firstpint=posp[j,1]
                            diff=firstr-firstp
                            aveB=(testDiff+posp[j,5])/2
                            testL=prBSpacing(aveB,gsK,(gsK+1))
                            #The spacing should be 4B for this set
                            if testL-tol <= diff <= testL+tol and lowIntTol <= (firstpint/firstrint) <= hiIntTol:
                                #This is a match we need to dump the lines into outArr
                                death=np.zeros([0,1])
                                match=zipSet(posr, i, posp,j,death,-1)
                                labArr=np.zeros([len(match),1])
                                labArr[:,0]=label
                                match=np.hstack([match,labArr])
                                outArr=np.vstack([outArr,match])
                                label+=1
            return outArr, inArr            
        else:
            posr,inr=findRBranch(inArr, gsK, wiggle,gsK)
            posp,inp=findPBranch(inArr, (gsK+1),wiggle,gsK)
            posq,inq=findQBranch(inArr, gsK,wiggle,gsK)
            
            for i in range(len(posr)):
                if posr[i,4] == gsK:
                    firstr=posr[i,0]
                    firstrint=posr[i,1]
                    #this was the matched 2B value
                    testDiff=posr[i,5]
                    for j in range(len(posp)):
                        if posp[j,4] == (gsK+1) and 0.975 <= (testDiff/posp[j,5]) <= 1.025:
                            firstp=posp[j,0]
                            firstpint=posp[j,1]
                            diff=firstr-firstp
                            aveB=(testDiff+posp[j,5])/2
                            testLr=prBSpacing(aveB,gsK,(gsK+1))
                            testLq=qBSpacing(aveB,gsK)
                            if testLr-tol <= diff <= testLr+tol and lowIntTol <= (firstpint/firstrint) <= hiIntTol:
                                #This is a match we need to dump the lines into outArr
                                for k in range(len(posq)):
                                    #firstQ should be 4B away for this set
                                    if posq[k,4] == gsK and testLq-(tol) <= (firstr-posq[k,0]) <= testLq+(tol) and lowIntTol/2 <= (posq[k,1]/firstrint) <= 2*hiIntTol:
                                        match=zipSet(posr, i, posp,j,posq,k)
                                        labArr=np.zeros([len(match),1])
                                        labArr[:,0]=label
                                        match=np.hstack([match,labArr])
                                        outArr=np.vstack([outArr,match])
                                        label+=1
            return outArr, inArr       
    if delK == 1:
        
            posr,inr=findRBranch(inArr, gsK, wiggle,gsK)
            posp,inp=findPBranch(inArr, (gsK+2),wiggle,gsK)
            posq,inq=findQBranch(inArr, gsK,wiggle,gsK) 
            for i in range(len(posr)):
                if posr[i,4] == gsK:
                    firstr=posr[i,0]
                    firstrint=posr[i,1]
                    #this was the matched 2B value
                    testDiff=posr[i,5]
                    for j in range(len(posp)):
                        if posp[j,4] == (gsK+2) and 0.975 <= (testDiff/posp[j,5]) <= 1.025:
                            firstp=posp[j,0]
                            firstpint=posp[j,1]
                            diff=firstr-firstp
                            aveB=(testDiff+posp[j,5])/2
                            testLr=prBSpacing(aveB,gsK,(gsK+2))
                            testLq=qBSpacing(aveB,gsK)
                            #The spacing should be 6B for this set
                            if testLr-tol <= diff <= testLr+tol and lowIntTol <= (firstpint/firstrint) <= hiIntTol:
                                #This is a match we need to dump the lines into outArr
                                for k in range(len(posq)):
                                    #firstQ should be 2B away for this set
                                    if posq[k,4] == gsK and testLq-(tol) <= (firstr-posq[k,0]) <= testLq+(tol) and lowIntTol/2 <= (posq[k,1]/firstrint) <= 2*hiIntTol:
                                        match=zipSet(posr, i, posp,j,posq,k)
                                        labArr=np.zeros([len(match),1])
                                        labArr[:,0]=label
                                        match=np.hstack([match,labArr])
                                        outArr=np.vstack([outArr,match])
                                        label+=1
            return outArr, inArr  
    if delK == -1:
        if 0 < gsK:
            posr,inr=findRBranch(inArr, gsK, wiggle,gsK)
            posp,inp=findPBranch(inArr, gsK, wiggle,gsK)
            posq,inq=findQBranch(inArr, gsK, wiggle,gsK)
            for i in range(len(posr)):
                if posr[i,4] == gsK:
                    firstr=posr[i,0]
                    firstrint=posr[i,1]
                    #this was the matched 2B value
                    testDiff=posr[i,5]
                    for j in range(len(posp)):
                        if posp[j,4] == gsK and 0.975 <= (testDiff/posp[j,5]) <= 1.025:
                            firstp=posp[j,0]
                            firstpint=posp[j,1]
                            diff=firstr-firstp
                            aveB=(testDiff+posp[j,5])/2
                            testLr=prBSpacing(aveB,gsK,gsK)
                            testLq=qBSpacing(aveB,gsK)
                            #The spacing should be 6B for this set
                            if testLr-tol <= diff <= testLr+tol and lowIntTol <= (firstpint/firstrint) <= hiIntTol:
                                #This is a match we need to dump the lines into outArr
                                for k in range(len(posq)):
                                    #firstQ should be 4B away for this set
                                    if posq[k,4] == gsK and testLq-(tol) <= (firstr-posq[k,0]) <= testLq+(tol) and lowIntTol/2 <= (posq[k,1]/firstrint) <= 2*hiIntTol:
                                        match=zipSet(posr, i, posp,j,posq,k)
                                        labArr=np.zeros([len(match),1])
                                        labArr[:,0]=label
                                        match=np.hstack([match,labArr])
                                        outArr=np.vstack([outArr,match])
                                        label+=1
            return outArr, inArr 
    else:
        print("You have violated selection rules\n")
        return 0

#Test point for stitchTog
#out, inp=stitchTog(data,0,0, diffTol)
#writeOut('test', out) 
        
#Ask the user what they want to look for
def getExclude():
    """Query to the user for values of delta K they do not want to search for
        
    Returns
    -------
    exList = the values of deltaK to not include
        
    """
    print('Enter any delK [1,0,-1] values that you do NOT want to search for: '+'\n')
    exList=[]
    print('Enter "done" if you are finished entering values'+'\n')
    while True:
        try:
            getV=input('Enter delK value: ')
            if getV == 'done':
                break
            add=int(getV)
            exList.append(add) 
        except ValueError:
            print('The value you entered in not valid'+'\n')
            continue   
    return exList
        
       
#loop through the spectrum and look for P/Q/R branches that match the pattern
#and fit the adjacent K sets  

def matchK(inArr, inmaxK,indelK,secTol):
    """Match the results of stitchTog agasint each other on the basis of the expect spacing between K stacks
        
    Parameters
    ----------
    inArr : array of the observed transitions
    
    inmaxK : the maximum value of K to test
    
    indelK : the delta K value to test
    
    secTol : the tolerance for the match defined by the user
    
    Returns
    -------
    outArr : the output matched array
    
    """
    #Need to generate an array for each gsK value and collect them in a master
    #array
    delK=int(indelK)
    outArr=np.zeros([0,10])
    masterArr=np.zeros([0,9])
    for i in range(0,inmaxK+1):
        if delK ==-1 and i ==0:
            #this is a selection rule thing
            continue
        else:
            tempArr,reArr=stitchTog(inArr,delK,i,secTol)
            masterArr=np.vstack([masterArr,tempArr])
            print('Found gsK= '+'%s'%i)
    print('Total Hits = '+'%s'%len(masterArr))
    #this is a test point, please don't uncomment
    #print(len(masterArr))
    #now need to loop through this array and pull out the sets that match the criteria
    usedLabels=[]
    for i in range(len(masterArr)):
        if masterArr[i,0]!= -1 and masterArr[i,8] not in usedLabels[:]:
            #this will be the array that collects the match
            #Get a set that we haven't tested before
            usedLabels.append(masterArr[i,8])
            #Now collect the set in dummmy array
            temp0=0
            ave0=0
            dummy=np.zeros([0,9])
            for l in range(len(masterArr)):
                if i <= l and masterArr[i,8] == masterArr[l,8]:
                    dummy=np.vstack([dummy,masterArr[l,:]])
                    if masterArr[l,6]!=0:
                        temp0+=1
                        ave0+=masterArr[l,5]
            #we now have an array "dummy" containing our set
            #need to collect our
            aveB0=(ave0/temp0)
            testedLabels=[]
            for j in range(len(masterArr)):
                colArr=np.zeros([0,9])
                if masterArr[j,0] != -1  and masterArr[j,8] not in testedLabels[:]:
                    if masterArr[i,8] != masterArr[j,8] and masterArr[i,7] != masterArr[j,7]:
                        colArr=np.vstack([colArr,dummy])
                        #now we can create a test set
                        dummy2=np.zeros([0,9])
                        testedLabels.append(masterArr[j,8])
                        temp=0
                        ave=0
                        for k in range(len(masterArr)):
                            if j <= k and masterArr[j,8] == masterArr[k,8]:
                                dummy2=np.vstack([dummy2,masterArr[k,:]])
                                if masterArr[k,6]!=0:
                                    ave+=masterArr[k,5]
                                    temp+=1
                        aveB1=(ave/temp)
                        colArr=np.vstack([colArr,dummy2])
                        if 0.975 <= (aveB1/aveB0) <= 1.025:
                            #now we have a array we can compare to dummy
                            #we need to test the r branch value only
                            #may add tests for P and Q branch later
                            #####
                            #Okay first we need to get common R lines for this 
                            #we need to figure out the J values of the first R lines
                            mult=dummy2[0,4]-dummy[0,4]
                            #figure out were the adjacent R line should be 
                            adjR=dummy[0,0]+(mult*aveB1)
                            testDiff=dummy2[0,0]-adjR
                            gsK0=dummy[0,7]
                            gsK1=dummy2[0,7]
                            mod=(gsK1-gsK0)
                            if indelK == 0:
                                #now solve for the tolerance
                                tol=0
                                tol=(testDiff/((2*gsK0*mod)+(mod**2)))
                                if -1*(inM/2) < tol < (inM/2):
                                    secKlist=[]
                                    secLablist=[]
                                    secKlist.append(gsK0)
                                    secKlist.append(gsK1)
                                    for x in range(len(masterArr)):
                                        if masterArr[x,0]!=-1 and masterArr[x,7]!=masterArr[j,7]:
                                            if masterArr[x,7] not in secKlist[:] and masterArr[x,8] not in secLablist[:]:
                                                dummy3=np.zeros([0,9])
                                                secLablist.append(masterArr[x,8])
                                                temp2=0
                                                ave2=0
                                                for y in range(len(masterArr)):
                                                    if x<=y and masterArr[y,8]==masterArr[x,8]:
                                                        dummy3=np.vstack([dummy3,masterArr[y,:]])
                                                        if masterArr[y,6]!=0:
                                                            ave2+=masterArr[y,5]
                                                            temp2+=1
                                                aveB2=(ave2/temp2)
                                                if 0.975<=(aveB2/aveB1)<=1.025:
                                                    #now we need to do the comparison
                                                    mult2=dummy3[0,4]-dummy2[0,4]
                                                    adjR2=dummy2[0,0]+(mult2*aveB2)
                                                    testD2=dummy3[0,0]-adjR2
                                                    gsK2=dummy3[0,7]
                                                    mod2=(gsK2-gsK1)
                                                    if 0.95*((2*tol*mod2*gsK1)+((mod2**2)*tol)) <= testD2 <= 1.05*((2*tol*mod2*gsK1)+((mod2**2)*tol)):
                                                        #this is a match yay
                                                        secKlist.append(gsK2)
                                                        colArr=np.vstack([colArr,dummy3])
                                    if 3 < len(secKlist):
                                        #another test point
                                        print(len(colArr))
                                        breaker=np.array([-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])
                                        tolin=np.zeros([len(colArr),1])
                                        tolin[:,0]=tol
                                        colArr=np.hstack([colArr,tolin])
                                        outArr=np.vstack([outArr,breaker])
                                        outArr=np.vstack([outArr, colArr])
                                    else:
                                        continue
                            if indelK == 1:
                                #now solve for the tolerance
                                tol=0
                                tol=((testDiff-(2*inM*mod))/((2*gsK1*mod)+(mod**2)))
                                if -1*(inM/2) < tol < (inM/2):
                                    secKlist=[]
                                    secLablist=[]
                                    secKlist.append(gsK0)
                                    secKlist.append(gsK1)
                                    for x in range(len(masterArr)):
                                        if masterArr[x,0]!=-1 and masterArr[x,7]!=masterArr[j,7]:
                                            if masterArr[x,7] not in secKlist[:] and masterArr[x,8] not in secLablist[:]:
                                                dummy3=np.zeros([0,9])
                                                secLablist.append(masterArr[x,8])
                                                temp2=0
                                                ave2=0
                                                for y in range(len(masterArr)):
                                                    if x<=y and masterArr[y,8]==masterArr[x,8]:
                                                        dummy3=np.vstack([dummy3,masterArr[y,:]])
                                                        if masterArr[y,6]!=0:
                                                            ave2+=masterArr[y,5]
                                                            temp2+=1
                                                aveB2=(ave2/temp2)
                                                if 0.975<=(aveB2/aveB1)<=1.025:
                                                    #now we need to do the comparison
                                                    mult2=dummy3[0,4]-dummy2[0,4]
                                                    adjR2=dummy2[0,0]+(mult2*aveB2)
                                                    testD2=dummy3[0,0]-adjR2
                                                    gsK2=dummy3[0,7]
                                                    mod2=(gsK2-gsK1)
                                                    if 0.95*((2*tol*mod2*gsK1)+((mod2**2)*tol)+(2*inM*mod2)) <= testD2 <= 1.05*((2*tol*mod2*gsK1)+((mod2**2)*tol)+(2*inM*mod2)):
                                                        #this is a match yay
                                                        secKlist.append(gsK2)
                                                        colArr=np.vstack([colArr,dummy3])
                                    if 3 < len(secKlist):
                                        #another test point
                                        print(len(colArr))
                                        breaker=np.array([-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])
                                        tolin=np.zeros([len(colArr),1])
                                        tolin[:,0]=tol
                                        colArr=np.hstack([colArr,tolin])
                                        outArr=np.vstack([outArr,breaker])
                                        outArr=np.vstack([outArr, colArr])
                                    else:
                                        continue
                                    
                            if indelK == -1:
                                #now solve for the tolerance
                                tol=0
                                tol=((testDiff+(2*inM*mod))/((2*gsK1*mod)+(mod**2)))
                                if -1*(inM/2) < tol < (inM/2):
                                    secKlist=[]
                                    secLablist=[]
                                    secKlist.append(gsK0)
                                    secKlist.append(gsK1)
                                    for x in range(len(masterArr)):
                                        if masterArr[x,0]!=-1 and masterArr[x,7]!=masterArr[j,7]:
                                            if masterArr[x,7] not in secKlist[:] and masterArr[x,8] not in secLablist[:]:
                                                dummy3=np.zeros([0,9])
                                                secLablist.append(masterArr[x,8])
                                                temp2=0
                                                ave2=0
                                                for y in range(len(masterArr)):
                                                    if x<=y and masterArr[y,8]==masterArr[x,8]:
                                                        dummy3=np.vstack([dummy3,masterArr[y,:]])
                                                        if masterArr[y,6]!=0:
                                                            ave2+=masterArr[y,5]
                                                            temp2+=1
                                                aveB2=(ave2/temp2)
                                                if 0.975<=(aveB2/aveB1)<=1.025:
                                                    #now we need to do the comparison
                                                    mult2=dummy3[0,4]-dummy2[0,4]
                                                    adjR2=dummy2[0,0]+(mult2*aveB2)
                                                    testD2=dummy3[0,0]-adjR2
                                                    gsK2=dummy3[0,7]
                                                    mod2=(gsK2-gsK1)
                                                    if 0.95*((2*tol*mod2*gsK1)+((mod2**2)*tol)-(2*inM*mod2)) <= testD2 <= 1.05*((2*tol*mod2*gsK1)+((mod2**2)*tol)-(2*inM*mod2)):
                                                        #this is a match yay
                                                        secKlist.append(gsK2)
                                                        colArr=np.vstack([colArr,dummy3])
                                    if 3 < len(secKlist):
                                        #another test point
                                        print(len(colArr))
                                        breaker=np.array([-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])
                                        tolin=np.zeros([len(colArr),1])
                                        tolin[:,0]=tol
                                        colArr=np.hstack([colArr,tolin])
                                        outArr=np.vstack([outArr,breaker])
                                        outArr=np.vstack([outArr, colArr])
                                    else:
                                        continue        
    return outArr
                                    
                    
        
                

       
 

################################################################################
#Ok now I think we can setup the main function
#What we want to do is to call stitchTog seqentially so that we search the entire
#spectrum, but we want to do this such that after a match is found we can delete those lines 
#and continue searching.

#We still have to chose where to start, ie what are the first inputs for stitchTog
#Normally I think I will start with delK=0, K"=0 eventhough that gives a large # of lines
#also remember that in order for a match to be found it must have 3 members in each allowed P/Q/R branch [working on getting around that]
def main():
    print("\nWelcome, just a reminder that the units you enter\nmust match the units of the input file's freq.\n")
    time.sleep(1)
    global upper2B 
    global lower2B
    global wiggle
    global lowIntTol
    global hiIntTol
    global diffTol
    global inM
    global label
    label=0
    sysname=input('What is the system you are looking at?  ')
    inArrname=input('What is the Line file name (including extension)?  ')
    rawinArr=np.loadtxt(inArrname)
    inArr=indexData(rawinArr)
    lowerB=input('Please enter lower B..  ')
    upperB=input('Please enter upper B..  ')
    lower2B=float(lowerB)*2
    upper2B=float(upperB)*2
    wiggle=float(input('Please enter tolerance for the 2B spacing (typically 50 Mhz/0.002cm-1)..  '))
    lowIntTol=float(input('Please enter lower intensity tolerance(typically=0.5)..  '))
    hiIntTol=float(input('Please enter higher intensity tolerance(typically=1.5)..  '))
    diffTol=float(input('Please enter tolerance for the P/Q/R spacing (typically 300-1000Mhz/0.01-0.033cm-1)..  '))
    version=input('Would you like to run advanced search? (y/n)   ')
    if version =='y':
        #this is the fancy version with the new function
        inM=float(input('Please enter the value for M..  '))
        maxK=int(input('Please enter the max ground state K value...  '))
        print("\nRunning...\n")
        
        for i in {1,-1,0}:
            print('Starting delK='+'%s'%i+'\n')
            outArr=matchK(inArr,maxK,i,diffTol)
            #this is a test point, please don't uncomment
            #print(outArr)
            if 0 < len(outArr):
                inArr=killLine(inArr,outArr)
                print('Found matches for delK='+'%s'%i+'\n')
                writeOut('delK_'+'%s'%i+'_fullmatch_',outArr,10)
        f=open('match_input_README'+'_'+'%s'%today+'.txt','w')
        f.write('These are the input parameters for pattern_match.py'+'\n')
        f.write('System =  '+'%s'%sysname+'\n')
        f.write('Filename =  '+'%s'%inArrname+'\n')
        f.write('Lower B value =  '+'%s'%lowerB+'\n')
        f.write('Upper B value =  '+'%s'%upperB+'\n')
        f.write('2B tolerance =  '+'%s'%wiggle+'\n')
        f.write('Lower intensity tolerance =  '+'%s'%lowIntTol+'\n')
        f.write('Upper intensity tolerance =  '+'%s'%hiIntTol+'\n')
        f.write('P/Q/R tolerance =  '+'%s'%diffTol+'\n')
        f.write('Max K value value searched =  '+'%s'%maxK+'\n')
        f.write('Value for M = '+'%s'%inM+'\n')
        f.close()
        f=open('Unassigned_lines'+'_'+'%s'%today+'.txt','w')
        for i in range(len(inArr)):
            for j in range(0,3):
                f.write('%s'%inArr[i,j]+'\t')
            f.write('\n')
        f.close()
        print("Have a nice day :)\n")
    #the else clause is the normal version
    else:
        startDelK=int(input('Please enter the starting del(K) value...  '))
        startgsK=int(input('Please enter the starting ground state K value...  '))
        maxK=int(input('Please enter the max ground state K value...  '))
        #Explicit print statement for input parameters
        f=open('match_input_README'+'_'+'%s'%today+'.txt','w')
        f.write('These are the input parameters for pattern_match.py'+'\n')
        f.write('System =  '+'%s'%sysname+'\n')
        f.write('Filename =  '+'%s'%inArrname+'\n')
        f.write('Lower B value =  '+'%s'%lowerB+'\n')
        f.write('Upper B value =  '+'%s'%upperB+'\n')
        f.write('2B tolerance =  '+'%s'%wiggle+'\n')
        f.write('Lower intensity tolerance =  '+'%s'%lowIntTol+'\n')
        f.write('Upper intensity tolerance =  '+'%s'%hiIntTol+'\n')
        f.write('P/Q/R tolerance =  '+'%s'%diffTol+'\n')
        f.write('Max K value value searched =  '+'%s'%maxK+'\n')
        f.close()
        print("\nRunning...\n")
        print("Starting delK="+'%s'%startDelK+', gsK='+'%s'%startgsK)
        ###First simple method
        firstout,firstin=stitchTog(inArr,startDelK,startgsK,diffTol)
        #now we want to remove lines that were assigned
        remainer=killLine(firstin,firstout)
        #build new array and writeOut and continue
        if 0 < len(firstout):
            print("Found matches for delK="+'%s'%startDelK+", gsK="+'%s'%startgsK)
            writeOut('delK_'+'%s'%startDelK+'_gsK_'+'%s'%startgsK,firstout,9)
            print("Lines remaining="+'%s'%len(remainer)+'\n')
        for i in {1,-1,0}:
            for j in range(maxK,-1,-1):
                #This is a selection rule thing
                if i == -1 and j == 0:
                    continue
                #You can change this to whatever set you start with
                if i ==startDelK and j==startgsK:
                    continue
                else:
                    print("Starting delK="+'%s'%i+', gsK='+'%s'%j)
                    newout,newin=stitchTog(remainer, i,j,diffTol)
                    #repeat the kill and write process
                    if 0 < len(newout):
                        print("Found matches for delK="+'%s'%i+", gsK="+'%s'%j)
                        writeOut('delK_'+'%s'%i+'_gsK_'+'%s'%j,newout,9)
                        remainer=killLine(newin,newout)
                        print("Lines remaining="+'%s'%len(remainer)+'\n')
        #after the process it would be a good idea to also print out the list of unassigned lines
        f=open('Unassigned_lines'+'_'+'%s'%today+'.txt','w')
        for i in range(len(remainer)):
            for j in range(0,3):
                f.write('%s'%remainer[i,j]+'\t')
            f.write('\n')
        f.close()
        print("Have a nice day :)\n")

if __name__ == "__main__":
    main()    
