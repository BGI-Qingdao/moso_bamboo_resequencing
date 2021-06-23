#Feb. 17, 2020: re-write the multi-allelic hx model
#Feb. 14, 2020: revise the weight in multi-allelic model
#Feb 5, 2020: Bug fix
#Nov 26, 2019: Modified output format
#Oct 25, 2019: Fixed a bug in B1
#Oct 20, 2019: software pass initial debugging & benchmarking
import numpy as np
import itertools as itr
from math import log,exp
from datetime import datetime
import re, sys, argparse


vlog = np.vectorize(log)
vexp = np.vectorize(exp)

def getAd(alpha):
    return(-log(alpha))


''' Read genome-wide substitution:polymorphism ratio
#format: N sub poly; no header
#0 for sub, 1 for poly
'''
def getPolyDen(spectfile):
    g=[]; G=[]
    with open(spectfile,'r') as spect:
        l=spect.next()
        l=l.strip().split('\t')
        N=int(l[0]); s=float(l[1]); p=float(l[2])
        print('Substitutions: %s ; polymorphisms: %s' %(s,p))
    try:
        assert s+p == 1
    except:
        print('s+p= %s + %s = %s != 1' % (s,p,s+p))
        sys.exit()
    g={(0,N):s,(1,N):p}
    G={(0,N):log(s), (1,N):log(p)}
    return(g,G,N)

'''Read unormalized neutral spect
##format: x n freq; no header
'''
def getLogSpect(spectfile,nosub,MAF):
    g={}; G={}; N=[]; checksum=0.
    with open(spectfile,'r') as spect:
        for l in spect:
            l=l.strip().split("\t")
            x=int(l[0]); n=int(l[1])
            f=float(l[2])
            g[(x,n)]=f
            G[(x,n)]=log(f)
            checksum += f
            N.append(n)
    N = list(set(N))
    print('Total prob is %f in the spectrum file.' %(checksum))
    if len(N) >= 2:
        print("Current implementation only supports uniform sample size. Please try again.")
        sys.exit()
    else:
        N=N[0]
    #Make sure there's no missing keys:
    fold=False
    for k in range(N+1):
        if (k,N) not in g:
            g[(k,N)] = 0
        elif MAF and k > N/2:
            fold = True
    #if fold;     
    if MAF and fold:
        fg = {}; fG = {}; checksum=0
        for i in range(N/2):
            if nosub and i==0:
                print('Skip substitutions')
                continue
            fg[(i,N)] = g[(i,N)] + g[(N-i,N)]
            fG[(i,N)] = log(fg[(i,N)])
            checksum += fg[(i,N)]
        if N%2 == 0:
            fg[(N/2,N)] = g[(N/2,N)]
            fG[(N/2,N)] = log(fg[(N/2,N)])
            checksum += fg[(N/2,N)]
        print('After folding, total probability is %s.'%(checksum))
        return(fg,fG,N)
    #if not folded:
    return(g,G,N) #spect, logspect, size

#global vairables
ComposSet = {}

#return collections of ordered 1D mk-long vectors 
#all vectors are of same length
def composition(mk,nk): 
    global ComposSet
    if mk < 1 or nk < 0:
        return []
    elif mk == 1:
        if (1,nk) not in ComposSet:
            ComposSet[(1,nk)] = [[nk]]
            return [[nk]]
        else:
            return ComposSet[(1,nk)]
    elif nk == 0:
        if (mk,0) not in ComposSet:
            Compos = []
            vec = [0]*mk
            ComposSet[(mk,0)] = [vec]
            return [vec]
        else:
            return ComposSet[(mk,0)]
    else:
        if (mk,nk) not in ComposSet:
            Compos = []
            for i in range(nk+1):
                tempComp = composition(mk-1,nk-i)
                for vec in tempComp:
                    vec2 = [i] + vec
                    if vec2 not in Compos:
                        Compos.append(vec2)
            ComposSet[(mk,nk)] = Compos
            return Compos
        else:
            return ComposSet[(mk,nk)]

#partition n into k non-empty classes in descending order
def partition(n, k, m=None):
    if m is None:
        m = n - k + 1 # maximum can be n - k + 1 since minimum is 1
    if k == 0:
        if n == 0:
            yield ()
        return
    for i in range(1, m + 1): # first could be from 1 to the maximum
        # the rest of the sum will be n - i among k - 1 elements with
        # maximum i
        for t in partition(n - i, k - 1, i):
            yield (i, ) + t

'''return multinomial look-up table for all possible combos'''
def getMNprobs(m,n,X):
	MNcompos = composition(m,n)
	logX = vlog(np.array(X))
	logIs = np.array([0]+[log(k) for k in range(1,n+1)])
	#such that log(k) = logIs[k]
	Probs = {}
	for K in MNcompos:
		assert len(K) == len(X)
		facKs = np.array([sum(logIs[:k+1]) for k in K])
		K = np.array(K)
		#print len(K),len(logX),len(facKs[:k])
		prob = sum(logIs) + K.dot(logX) - facKs.dot(np.ones(len(K)))
		Probs[tuple(K)]=exp(prob)
	return(Probs)

def choose(n,k):
    if n <= k :
        return 1
    elif n <=0 or k <= 0:
        return 1
    logIs = [log(i) for i in range(1,n+1)]
    logN = sum(logIs) - sum(logIs[:k]) - sum(logIs[:(n-k)])
    return(exp(logN))
    
'''
m --> # of classes
n --> sample size
k --> # of derived alleles
X --> m-element vector for equilibrium frequency
return unnormalized probability
'''
def getKprob(m,n,X,k,nofreq=False):
    xProbs = {} ; Xpool={}
    for mi in range(1,m):
        for xi in itr.combinations(range(len(X)),r=mi):
            x = [X[i] for i in xi]
            xa = [X[i] for i in range(len(X)) if i not in xi]
            xt = np.array(x+[sum(xa)])
            xt = tuple(xt/sum(xt))
            xProbs[xt] = getMNprobs(len(x)+1,n,xt)
    prob_k=0 
    for mi in range(1,m):
        if mi > k:
            continue
        prob_j=0
        #for x in itr.combinations(X,r=mi):
        for xi in itr.combinations(range(len(X)),r=mi):
            x = [X[i] for i in xi]
            xa = [X[i] for i in range(len(X)) if i not in xi]
            xt = np.array(x+[sum(xa)])
            xt = tuple(xt/sum(xt))
            for kd in partition(k,mi):
                kt = list(kd) + [n-k]
                prob_j += xProbs[xt][tuple(kt)]
        prob_k += prob_j/choose(m,mi)
    return(prob_k)

'''
m --> # of classes
n --> sample size
k --> # of derived alleles
X --> m-element vector for equilibrium frequency
return full unnormalized probability
'''
def getFullKprobs(m,n,X):
    xProbs = {} 
    Xpool={}
    for mi in range(1,m):
        for xi in itr.combinations(range(len(X)),r=mi):
            x = [X[i] for i in xi]
            xa = [X[i] for i in range(len(X)) if i not in xi]
            xt = np.array(x+[sum(xa)])
            xt = tuple(xt/sum(xt))
            xProbs[xt] = getMNprobs(len(x)+1,n,xt)
    fullProbs={}
    for k in range(1,n+1):
        prob_k=0 
        for mi in range(1,m):
            if mi > k:
                continue
            prob_j=0
            #for x in itr.combinations(X,r=mi):
            for xi in itr.combinations(range(len(X)),r=mi):
                x = [X[i] for i in xi]
                xa = [X[i] for i in range(len(X)) if i not in xi]
                xt = np.array(x+[sum(xa)])
                xt = tuple(xt/sum(xt))
                for kd in partition(k,mi):
                    kt = list(kd) + [n-k]
                    prob_j += xProbs[xt][tuple(kt)]
            prob_k += prob_j/choose(m,mi)
        fullProbs[k] = prob_k
    return(fullProbs)

#just for B2 and B0; B1 normalization integrated in pipeline
def normalize(n,fullProbs,MAF=False, nosub=False):
	#base = sum(fullProbs.values())
	normProb = {}
	if nosub:
		del fullProbs[n]
	base = sum(fullProbs.values())
	if MAF:
		if not nosub:
			normProb[0] = fullProbs[n] / base
		for k in range(1,n/2):
			normProb[k] = ( fullProbs[k] + fullProbs[n-k] ) / base
		if n % 2 == 0:
			normProb[n/2] = fullProbs[n/2] / base
	else:
		for k in range(1,n):
			normProb[k] = fullProbs[k] / base
		if not nosub:
			normProb[n] = fullProbs[n] / base
	return(normProb)


#parse m and X
def parseXnM(m, argX):
    Xs = argX.split(',')
    if len(Xs) == 1:
        try:
            x=float(Xs)
            assert m == 2
            x=max(1-x, x)
            X=tuple(1-x,x)
            parsedX = [X]
        except(AssertionError):
            print('Please make sure vector X (via --fixX) contain the same number of values as m (via -m), the number of balanced alleles..')
            sys.exit()
        except: #X not entirely numerical
            s = Xs.replace('(','').replace(')','')
            s = tuple([float(_) for _ in s.split('|')])
            try:
                assert len(s) == m
            except(AssertionError):
                print('Please make sure vector X (via --fixX) contain the same number of values as m (via -m), the number of balanced alleles..')
                sys.exit()
            try:
                assert sum(s) == 1.
            except:
                print('Please make sure all values in vector X (via --fixX) sum up to 1.')
                sys.exit()
            parsedX = [s]
    else: #more than one value
        parsedX = []
        for X in Xs:
            try:
                x=float(Xs)
                assert m == 2
                x=max(1-x, x)
                X=tuple(1-x,x)
                parsedX.append(X)
            except(AssertionError):
                print('Please make sure vector X (via --fixX) contain the same number of values as m (via -m), the number of balanced alleles.')
                sys.exit()
            except: #X not entirely numerical
                s = Xs.replace('(','').replace(')','')
                s = tuple([float(_) for _ in s.split('|')])
                try:
                    assert len(s) == m
                except(AssertionError):
                    print('Please make sure the number of values provided via --fixX matches the presumed number of balanced alleles (-m).')
                    sys.exit()
                parsedX.append(s)
    return(parsedX,m)

'''Generate grids of parameters to optimize over'''
def getGrids(m, X, seqA, listA):
    #get X
    if X:
        Xgrid,m = parseXnM(m,X)
    else:
        #xGrid = [.05*i for i in range(1,11)]
        Xgrid = partition(20,m)
        Xgrid = [ tuple([0.05*i for i in parts]) for parts in Xgrid ]
    #get A
    if listA:
        AGrid = [float(x) for x in listA.split(',')]
    elif seqA:
        Amin,Amax,Astep = [float(x) for x in seqA.split(',')]
        n=(Amax-Amin)/Astep
        AGrid=[Amin+Atep*i for i in range(n+1)]
    else:
        AGrid = [100*i for i in range(1,12)] + [200*i for i in range(6,13)] + [500*i for i in range(5,10)] + [1000*i for i in range(5,11)] + [1e6,1e9]
    #define the default list of alpha values:
    aGrid = [0,1e-8,1e-7,1e-6,1e-5]+[m*1e-4 for m in range(1,10001)]
    #generate sorted AdGrid based on the dense grid of alpha:
    AdGrid = [-log(1e-32)]+[ getAd(a) for a in aGrid[1:] ] 
    return(Xgrid,AGrid,AdGrid,aGrid)
        

'''Generate look-up tables for g(k) and f(k,d; x,A) distributions
objects are logSpec[(k,N)] and Fx[X][a][k] respectively'''
def initialize(m,spectfile,xGrid,aGrid,MAF,nofreq,nosub):
    #Read spect and get sample size
    if nofreq:
        print('Only using polymorphism/substitution calls. 0 for substitution, 1 for polymorphism.')
        spec,logSpec,N = getPolyDen(spectfile)
        #print(spec)
        print('Sample size: %s'%(N))
        Fx={}
        for X in xGrid:
            p0 = getKprob(m,N,X,0,nofreq=True) #p0==pn
            pn = getKprob(m,N,X,N,nofreq=True)
            sub = pn/(1-p0); poly = 1-sub
            Hx = [sub, poly]
            Fx[X] = {}
            for a in aGrid:
                Fx[X][a] = {}
                for k in [0,1]: #0 for sub, 1 for poly
                    gx = spec[(k,N)]
                    fx = a*Hx[k] + (1-a)*gx
                    Fx[X][a][k] = log(fx)
    else:
        spec,logSpec,N = getLogSpect(spectfile,nosub,MAF)
        #print(spec)
        print('Sample size: %s'%(N))
        #get the range of k
        if MAF:
            print('Using minor allele frequencies...')
            kpool = list(range(N/2+1))
            if nosub:
                print('...without substitutions')
                kpool = list(range(1,N/2+1))
        else:
            print('Using derived allele frequencies...')
            kpool = list(range(1,N+1))
            if nosub:
                print('...without substitutions')
                kpool = list(range(1,N))
        #Generate unnormalized h(x)
        Fx = {}
        for X in xGrid:
            Hx = getFullKprobs(m,N,X)
            Hx = normalize(N,Hx,MAF=MAF, nosub=nosub)
            Fx[X] = {}
            #print(aGrid)
            for a in aGrid:
                #print(a)
                Fx[X][a] = {}
                for k in kpool:
                    gx = spec[(k,N)]
                    fx = a*Hx[k] + (1-a)*gx
                    try:
                        Fx[X][a][k] = log(fx)
                    except(ValueError):
                        print('Skipping (k=%s,N=%s) in full model for X=(%s), alpha:%g , fx: %g, hx: %g, gx:%g'%(k, N, ','.join([str(x) for x in X]),a, fx, Hx[k], gx))
                        continue
    return(spec,logSpec,Fx,N)


''' Return parsed data
# all sitePos are genetic positions (by cM)
# default value of Rrate is 1e-6 cM/site
# Input format: physPos, genPos, x, n
'''
def readInput(infile,nosub,nofreq,MAF,phys=False,Rrate=1e-6):
    phys_pos = []
    Ks=[]; Ns = []; pos_list=[]; numSites=0 
    postype=1-int(phys) #index; 0 is physical position, 1 is genetic position
    with open(infile,'r') as sites:
        l=sites.next()#skipping header
        if nosub:
            #make sure to only read frequencies
            for l in sites:
                l=l.strip().split('\t')
                if float(l[2])/float(l[3]) not in [0.,1.]:
                    numSites += 1
                    physPos,k,n = [float(l[0]),int(l[2]),int(l[3])]
                    #genPos = float(l[1])
                    physPos = int(float(l[0]))
                    sitepos = float(l[postype])*phys*Rrate + float(l[postype])*(1-int(phys))
                    k = MAF*min(k,n-k) + (1-MAF)*k
                    Ks.append(k); Ns.append(n); pos_list.append(sitepos)
                    phys_pos.append(physPos)
        elif nofreq: # Take all none N cases as 1
            for l in sites:
                l=l.strip().split('\t')
                numSites += 1
                physPos,k,n = [float(l[0]),int(l[2]),int(l[3])]
                #genPos = float(l[1])
                physPos = int(float(l[0]))
                sitepos = float(l[postype])*phys*Rrate + float(l[postype])*(1-int(phys))
                k = (k!=n) #k=1 if not sub, k=0 if sub
                Ks.append(k); Ns.append(n); pos_list.append(sitepos)
                phys_pos.append(physPos)
        else:
            for l in sites:
                l=l.strip().split('\t')
                numSites += 1
                physPos,k,n = [float(l[0]),int(l[2]),int(l[3])]
                #genPos = float(l[1])
                physPos = int(float(l[0]))
                sitepos = float(l[postype])*phys*Rrate + float(l[postype])*(1-int(phys))
                k = MAF*min(k,n-k) + (1-MAF)*k
                Ks.append(k); Ns.append(n); pos_list.append(sitepos)
                phys_pos.append(physPos)
    return(phys_pos,pos_list,Ks,Ns,numSites)


''' Return log likelihood for the given window
# testSite is genetic position (in recomb unit)
# logP_neut = logSpect[(k,N)]; logP_alt = Fx[x][a][k]'''
def calcBaller(pos_list,Ks,Ns, testSite, testSite_i, logSpect, Fx, xGrid, AGrid, AdGrid, aGrid,MAF):
    #note that the test site is included in Ks
    #Optimize over xGrid and AGrid
    L = len(AdGrid) ; numSites = len(pos_list)
    #testSite = pos_list[testSite_i]
    Tmax = [-100,0,0] #LR, x, A
    # Go through all A values
    for A in AGrid:
        # i indexes positions. pos_list ascending
        # c indexes Ad value in the pre-computed descendingly sorted list
        i_index=[] ; c_index=[]
        # going through sites from the center to outward
        i = max(testSite_i-1,0) ; c = len(AdGrid)-1
        # Leftward:
        while i >= 0 and c > 0 and testSite_i != 0:
            pos = pos_list[i]
            dist = testSite - pos
            Ad = float(A)*dist
            while AdGrid[c] < Ad and c >= 0:
                c -= 1
            #now AdGrid[c] >= Ad or c=0
            if c >=0:
                i_index.append(i)
                c_index.append(c)
            i -= 1
        # now i==0. Starting rightward from center
        #i=min(testSite_i+1,pos_list[-1]) 
        i = min(testSite_i+1,numSites) ; c = len(AdGrid)-1
        while i < len(pos_list) and c > 0:
            pos = pos_list[i]
            dist = pos - testSite
            Ad = float(A)*dist
            while AdGrid[c] < Ad and c >=0:
                c -= 1
            if c >= 0:
                i_index.append(i)
                c_index.append(c)
            i += 1
        # if noCenter, consider testSite_i too
        if testSite != pos_list[testSite_i]:
            pos = pos_list[testSite_i]
            dist = abs(pos - testSite)
            Ad = float(A)*dist
            while AdGrid[c] < Ad and c >=0:
                c -= 1
            if c>=0:
                i_index.append(testSite_i)
                c_index.append(c)
        #go through the grid of x
        #note that this grid only spans 0-0.5 
        for x in xGrid:
            La=0; L0=0
            #go through all the sites
            for j in xrange(len(i_index)):
                i = i_index[j] #index in the pos_list
                c = c_index[j]
                alpha = aGrid[c]
                k = Ks[i] ; N = Ns[i]
                La += Fx[x][alpha][k]
                L0 += logSpect[(k,N)]
            #get LR
            T = 2*(La - L0)
            if T >=Tmax[0]:
                Tmax = [T,x,A]
                winSize = len(i_index)
    return(Tmax,winSize)

'''Perform the scan with fixed window size of w (bp), with step size s (bp).'''
def scan_fixSize_noCenter(xGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec,Fx,N,Rrate,w,s,MAF=False):
    #w = float(r) 
    print("writing output to %s" % (outfile))
    with open(outfile,'w') as scores:
        scores.write('physPos\tgenPos\tLR\txhat\tAhat\tnumSites\n')
        start = int( floor(2*float(phys_pos[0])/w)*(w/2) )
        end = start + s ; midpos = start + s/2
        start_i=0; end_i=0; pos_i=0
        while midpos <= phys_pos[-1]:
            #define the window
            while phys_pos[start_i] < start:
                start_i +=1 
            while phys_pos[pos_i] < midpos:
                pos_i += 1
            while end_i < numSites:
                if phys_pos[end_i] < end:
                    end_i += 1
                else:
                    break
            if start_i >= end_i:
                scores.write('%g\t%s\t0\tNA\tNA\t0\n' % (midpos,midpos*Rrate))
                start+=s; midpos+=s; end+=s
            else:
                #Window = Ks[start_i:end_i+1]
                WinKs = Ks[start_i:pos_i] + [0] + Ks[pos_i:end_i+1]
                WinNs = Ns[start_i:pos_i] + Ns[pos_i] + Ns[pos_i:end_i+1]
                winPosList = pos_list[start_i:pos_i] + [midpos*Rrate] + pos_list[pos_i:end_i+1]
                #calcBaller args: pos_list,Ks,Ns, testSite, testSite_i, logSpect, biFx, xGrid, AGrid, AdGrid, aGrid,MAF
                #pos_list and phys_pos has matching indice
                Tmax,winSites = calcBaller(winPosList,WinKs, WinNs, midpos*Rrate, pos_i-start_i, logSpec, Fx, xGrid, AGrid, AdGrid, aGrid, MAF)
                scores.write('%g\t%s\t%s\t%s\t%g\t%d\n' % (midpos, midpos*Rrate, Tmax[0], ','.join(['%g'%(_) for _ in Tmax[1]]), Tmax[2], winSites))#
                #read in, take next step
                start+=s; midpos+=s; end+=s 
    scores.close()
    return(0)

'''Perform the scan with fixed window size of w (bp), centered on every s informative sites'''
def scan_fixSize_siteCenter(xGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec,Fx,N,Rrate,w,s,MAF=False):
    #w = float(r) 
    print("writing output to %s" % (outfile))
    with open(outfile,'w') as scores:
        scores.write('physPos\tgenPos\tLR\txhat\tAhat\tnumSites\n')
        i=0; 
        start_i=0; end_i=0
        while i < numSites:
            testSite=phys_pos[i]
            start =  max(0, testSite-r/2); end = min(testSite+r/2,phys_pos[-1])
            while phys_pos[start_i] < start:
                start_i += 1
            while end_i < numSites:
                if phys_pos[end_i] < end:
                    end_i += 1
                else:
                    break
            assert end_i >= start_i
            end_i = min(end_i, numSites-1)
            #Window = Ks[start_i:end_i+1]
            Tmax,winSites = calcBaller(pos_list[start_i:end_i+1], Ks[start_i:end_i+1], Ns[start_i:end_i+1], pos_list[i], i-start_i , logSpec, Fx, xGrid, AGrid, AdGrid, aGrid, MAF)
            scores.write('%s\t%s\t%s\t%s\t%g\t%d\n'% (phys_pos[i], pos_list[i], Tmax[0], ','.join(['%g'%(_) for _ in Tmax[1]]), Tmax[2], winSites))#
            i+=int(s)
    scores.close()
    return(0)

'''Perform the scan with site-based window, with s sites on either side'''
def scan_siteBased(xGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec,Fx,N,phys,Rrate,r,s=1,MAF=False):
    print("writing output to %s" % (outfile))
    with open(outfile,'w') as scores:
        scores.write('physPos\tgenPos\tLR\txhat\tAhat\tnumSites\n')
        i=0
        while i < numSites:
            testSite = pos_list[i]
            start_i = max(0, i-r) ; end_i = min(numSites,i+r+1)
            Tmax,winSites = calcBaller(pos_list[start_i:end_i], Ks[start_i:end_i], Ns[start_i:end_i], testSite, i-start_i , logSpec, Fx, xGrid, AGrid, AdGrid, aGrid, MAF)
            scores.write('%s\t%s\t%s\t%s\t%g\t%d\n' % (phys_pos[i], pos_list[i], Tmax[0], ','.join(['%g'%(_) for _ in Tmax[1]]), Tmax[2], winSites))
            i+=s
    scores.close()
    return(0)

def scan_alpha(xGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec,logSpec,Fx,N,s=1,MAF=False):
    print("writing output to %s" % (outfile))
    with open(outfile,'w') as scores:
        scores.write('physPos\tgenPos\tLR\txhat\tAhat\tnumSites\n')
        i=0
        while i < numSites:
            testSite = pos_list[int(i)]
            Tmax,winSites = calcBaller(pos_list, Ks, Ns, testSite, int(i), logSpec, Fx, xGrid, AGrid, AdGrid, aGrid, MAF)
            scores.write('%s\t%s\t%s\t%s\t%g\t%d\n' % (phys_pos[i], pos_list[i], Tmax[0], ','.join(['%g'%(_) for _ in Tmax[1]]), Tmax[2], winSites))
            i+=int(s)
    scores.close()
    return(0)


def scan(xGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec,Fx,N,size=False,r=0,s=1,phys=False,nofreq=False,MAF=False,noCenter=False,Rrate=1e-6):
    if size:
        print('You\'ve chosen to fix the scan window size.')
        if r==0:
            print('Please set a window width in bp with \"-w\" or \"--window\" command.')
            sys.exit()
        if not phys:
            print('Please make sure to use physical positions as coordinates if fixed-length windows are chosen. Scan will continue with physical positions.')
            phys = True
            #sys.exit()
        if noCenter:
            w = float(r) 
            print('Computing LR on %.3f kb windows on every %s bp. Using physical positions by default.' % (w/1e3, s))
            scan_fixSize_noCenter(xGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec,logSpec,Fx,N,Rrate,w,s,MAF)
        else: #site-centered
            w = float(r) 
            print('Computing LR on %.2f kb windows on every %g informative sites. Using physical positions by default.' % (w/1e3,s))
            scan_fixSize_siteCenter(xGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec,Fx,N,Rrate,w,s,MAF)
    # When fixSize == False, and radius (-r) provided. Scan with fixed number of sites
    elif r != 0: 
        print('Computing LR on every %s site/s, with %s informative sites on either side.' % (s, r))
        scan_siteBased(xGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec,Fx,N,phys,Rrate,r,s,MAF)
    #window size not given 
    #then use all data (but the test site)
    else:
        print('Computing LR on every %s site/s, using all the data with alpha >= 1e-8.' % (s))
        scan_alpha(xGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec,logSpec,Fx,N,s,MAF)
    print(str(datetime.now())+'. Scan finished.')


#give the config file given the concatenated input
def getConfig(infile,configfile):
    Config={}; numSites=0# N: [s,p]
    with open(infile,'r') as sites:
        l=next(sites)#skip the header by default
        for l in sites:
            x,n = [int(x) for x in l.strip().split('\t')[2:] ]
            if x==0:
                print('Please make sure the input has derived allele frequency. Sites with 0 observed allele count (k=0) will be ignored.\n')
                continue
            if n not in Config:
                Config[n] = [0,0]
            Config[n][0] += int(x==n)
            Config[n][1] += 1-int(x==n)
            numSites+=1
    sizes = sorted(Config.keys())
    with open(configfile,'w') as config:
        for N in sizes:
            config.write('%s\t%s\t%s\n' % ( N, Config[N][0]/float(numSites) , Config[N][1]/float(numSites) ))
    sites.close(); config.close()
    print('Done')

#give spectrum file given the concatenated input
def getSpect(infile,spectfile,nosub=False,MAF=False):
    Spect={}; numSites=0
    if nosub:
        print('Generating spectrum for all polymorphic sites. Substitutions (x=n or x=0) won\'t be considered.')
        with open(infile,'r') as sites:
            l=next(sites)
            for l in sites:
                (x,n)=[ int(i) for i in l.strip().split('\t')[2:] ]
                if MAF:
                    x = min(x, n-x)
                if x == 0 or x == n:
                    continue
                if (x,n) in Spect:
                    Spect[(x,n)] += 1
                else:
                    Spect[(x,n)] = 1
                numSites+=1
    else:
        with open(infile,'r') as sites:
            l=next(sites)
            for l in sites:
                (x,n)=[ int(i) for i in l.strip().split('\t')[2:] ]
                if MAF:
                    x = min(x, n-x)
                elif x==0:
                    print('Please make sure the input has derived allele frequency. Sites with 0 observed allele count (k=0) should not be included.\n')
                    sys.exit()

                if (x,n) in Spect:
                    Spect[(x,n)] += 1
                else:
                    Spect[(x,n)] = 1
                numSites+=1
    #write out
    pairs = sorted(Spect.keys())
    with open(spectfile,'w') as spec:
        for x,n in pairs:
            spec.write('%s\t%s\t%s\n' % (x,n,float(Spect[(x,n)])/float(numSites) ))
    sites.close(); spec.close()
    print('Done.')

#main function to scan through the input
#only work with minor allele frequency
def main(): 
    #parsing arguments
    parser = argparse.ArgumentParser()#usage='python %prog -i <input file> -o <output file> --spect <spect/config file> [--help] [--nofreq] [--nosub] [--MAF] [--getSpect] [--getConfig] [--fixSize] [--physPos] [--rec <recomb rate>] [-w <window size>] [--noCenter] [-s <step size>] [--rangeA <min,max,step>] [--listA <A1,A2,..,Ak>] [--fixX <x>]'
    parser.add_argument('-i','--input', dest='infile', help = 'Path and name of your input file.\n',required=True)
    parser.add_argument('--spect',dest='spectfile', help = 'Path and name of the allele frequency spectrum file or configuration file.\n',required=True)
    parser.add_argument('-o','--output', dest='outfile', help = 'Path and name of your output file.\n')
    parser.add_argument('-m', dest='m', default=2, help = 'Presumed number of balanced alleles. Default is 2.\n')
    parser.add_argument('--getSpect', dest='getSpec', action='store_true', default=False, help='Option to generate frequency spectrum file from the concatenated input file. Use \"-i\" and \"--spect\" commands to provide names and paths to input and output files, respectively. Indicate the input type with \"--MAF\" and/or \"--nosub\".\n')
    parser.add_argument('--getConfig', dest='getConfig', action='store_true', default=False, help='Option to generate configuration file from the concatenated input file. Use \"-i\" and \"--spect\" commands to provide names and paths to input and output files, respectively.\n')

    parser.add_argument('--nofreq', dest='nofreq', action='store_true', default=False, help = 'Option to ignore allele frequency information. All polymorphic sites will be considered as equivalent.')
    parser.add_argument('--nosub', dest='nosub', action='store_true', default=False, help = 'Option to not include substitution in input data.')
    parser.add_argument('--MAF', dest='MAF', action='store_true', default=False, help = 'Option to use minor allele frequency, instead of polarized allele frequency. The latter is default.')
    parser.add_argument('--physPos', action='store_true', dest = 'phys', default = False, help = 'Option to use physical positions instead of genetic positions (in cM). Default is using genetic positions.\n')
    parser.add_argument('--rec', dest='Rrate', type = float, default = 1e-6 , help='The uniform recombination rate in cM/nt. Default value is 1e-6 cM/nt. Only useful when choose to use physical positions as coordinates.\n\n')

    parser.add_argument('--fixSize', action='store_true', dest = 'size', default = False, help = 'Option to fix the size of scanning windows. When true, provide the length of window in neucleotide (nt) with \"-w\" or \"--window\" command.\n')
    parser.add_argument('-w','--window', dest='r', type = int, default=0, help='Number of sites flanking the test locus on either side. When choose to fix window size (\"--fixSize\"), input the length of window in bp.\n')
    parser.add_argument('--noCenter', action='store_true', dest='noCenter', default=False, help = 'Option to have the scanning windows not centered on informative sites. Require that the window size (\"-w\") in physical positions (\"--physPos\") is provided. Default is True.\n')
    parser.add_argument('-s','--step', dest='step', type = float, default=1, help='Step size in bp (when using \"--noCenter\") or the number of informative sites. Default value is one site or one nucleotide.\n\n')

    parser.add_argument('--fixX', dest='x', help='Option to fix the presumed equilibrium frequency. When m>2, use the format \"(x1|x2|...|xm),(x1|x2|...|xm),...\" with each set of x1, x2,..., xm presented in descending order and summed up to one. \n')
    parser.add_argument('--rangeA', dest='seqA', help='Range of the values of the parameter A to optimize over. Format should follow <Amin>,<Amax>,<Astep> with no space around commas.\n')
    parser.add_argument('--listA', dest='listA', help='Manually provide a list of A values to optimize over. Please separate the values with comma, no space.\n')

    
    #if len(sys.argv[1:]) == 0:
    #    parser.print_help()
    #    sys.exit()

    opt = parser.parse_args()

    if opt.getSpec:
        print('You\'ve chosen to generate site frequency spectrum...')
        print('Concatenated input: %s \nSpectrum file: %s' % (opt.infile, opt.spectfile ))
        getSpect(opt.infile, opt.spectfile, opt.nosub, opt.MAF)
        sys.exit()
    elif opt.getConfig:
        print('You\'ve chosen to generate the substitution-polymorphism configuration...')
        print('Concatenated input: %s \nConfiguration file: %s' % (opt.infile, opt.spectfile ))
        getConfig(opt.infile, opt.spectfile)
        sys.exit()

    #generate the grids to optimize over/with
    xGrid,AGrid,AdGrid,aGrid = getGrids(int(opt.m), opt.x, opt.seqA, opt.listA)
    #print xGrid
    print('\nOptimizing over x= '+'; '.join( ["("+','.join(['%g'%(_) for _ in x])+")" for x in xGrid] ) +'.\n' ) 
    print('Optimizing over A= '+', '.join([str(a) for a in AGrid])+'.\n' ) 

    #initialization for a grid of x
    print('\n%s. Initializing...'%(datetime.now()))
    (Spec, logSpec, Fx,N) = initialize(int(opt.m),opt.spectfile,xGrid,aGrid,opt.MAF,opt.nofreq,opt.nosub)
    #start reading data
    print("\n%s. Reading input file: %s" % (datetime.now(),opt.infile) )
    (phys_pos,pos_list,Ks,Ns,numSites) = readInput(opt.infile,opt.nosub,opt.nofreq,opt.MAF,opt.phys,opt.Rrate)

    #finished reading file, start writing output
    print("\n%s. Start computing likelihood ratios..." %(datetime.now()))
    scan(xGrid,AGrid,AdGrid,aGrid,opt.outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec, Fx,N, opt.size,opt.r,opt.step,opt.phys,opt.nofreq,opt.MAF,opt.noCenter,opt.Rrate)#, ProbsStar

    #pipeline finished
    print('\n%s. Pipeline finished.'%(datetime.now()))



if __name__ == '__main__':
    main()
