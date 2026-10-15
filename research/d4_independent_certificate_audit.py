"""Second implementation of the d=4 three-copy certificate audit.

Imports no construction or previous verification program. Input archives are
read only. Literal permutations, partial traces, source generators, all Gram
maps and the final rational identity are rebuilt. Not an external peer review.
"""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product, combinations_with_replacement
import hashlib
import json
import math
import multiprocessing as mp
from pathlib import Path
import time
import numpy as np
from scipy import sparse
import sympy as sp
from flint import fmpq, fmpq_mat, nmod_mat

ROOT = Path(__file__).resolve().parent / 'matched_d4'
RUN = ROOT / 'full_gram_crossed'
CERT = RUN / 'rational_certificate'
OUT = RUN / 'audit_2026-10-14'
P = tuple(permutations(range(4)))
P3 = tuple(permutations(range(3)))
PI = {p:i for i,p in enumerate(P)}
MUL = np.array([[PI[tuple(p[q[i]] for i in range(4))] for q in P] for p in P])
INV = [PI[tuple(p.index(i) for i in range(4))] for p in P]
AUX = (0,3,4,7,8,9,11,12,13,15,16,18,20,23)
FAMILIES = ('31','22','13')

def read(path):
    return json.loads(path.read_text())

def rational(rows):
    return sp.Matrix([[sp.Rational(v) for v in row] for row in rows])

@lru_cache(None)
def literal(d, n=4, cut=()):
    # Start with the identity tensor. Transpose its OUTPUT axes by p^{-1}.
    ident = np.eye(d**n,dtype=np.int64).reshape((d,)*(2*n))
    result=[]
    for p in permutations(range(n)):
        order=[p.index(i) for i in range(n)]+list(range(n,2*n))
        a=ident.transpose(order)
        order=list(range(2*n))
        for i in cut:order[i],order[n+i]=order[n+i],order[i]
        result.append(a.transpose(order).reshape(d**n,d**n))
    return np.array(result)

def trace_slot(a,d,slot,order):
    b=np.trace(a.reshape((d,)*8),axis1=slot,axis2=4+slot)
    remaining=[j for j in range(4) if j!=slot]
    perm=[remaining.index(j) for j in order]
    return b.transpose(perm+[j+3 for j in perm]).reshape(d**3,d**3)

def insert_identity(a,d,slot):
    # Preserve the original order of the other three slots.
    b=np.kron(np.eye(d,dtype=np.int64),a).reshape((d,)*8)
    order=list(range(1,4));order.insert(slot,0)
    return b.transpose(order+[j+4 for j in order]).reshape(d**4,d**4)

def quotient(diagrams,keep):
    vectors=sp.Matrix(diagrams.reshape(len(diagrams),-1).T.tolist())
    basis=vectors[:,list(keep)]
    pivots=list(basis.T.rref()[1])
    left=basis[pivots,:].inv()
    q=left*vectors[pivots,:]
    assert basis*q==vectors
    return q,pivots,left

def local_data():
    aux,piv4,left4=quotient(literal(2),AUX)
    _,pivp,leftp=quotient(literal(4,3),range(6))
    aux3,piv3,left3=quotient(literal(2,3),range(5))
    def coordinates(a,piv,left):
        return list(left*sp.Matrix(a.ravel()[piv].tolist()))
    conditional={};marginals={}
    for slot in (0,3):
        order=[j for j in range(4) if j!=slot]
        phys=[coordinates(trace_slot(a,4,slot,order),pivp,leftp) for a in literal(4)]
        aq=[coordinates(2*a-insert_identity(trace_slot(a,2,slot,order),2,slot),piv4,left4)
            for a in literal(2)]
        conditional[slot]=(phys,[[v/2 for v in row] for row in aq])
    for name,slot,order in [('aab',2,[0,1,3]),('abb22',1,[0,2,3]),('abb13',1,[3,0,2])]:
        phys=[coordinates(trace_slot(a,4,slot,order),pivp,leftp) for a in literal(4)]
        aq=[coordinates(trace_slot(a,2,slot,order),piv3,left3) for a in literal(2)]
        marginals[name]=(phys,aq)
    return aux,conditional,marginals

def group(f):
    return [i for i,p in enumerate(P) if (p[3]==3 if f!='22' else set(p[:2])=={0,1})]

def generator(word,f):
    # Enumerate BOTH physical permutations and the two sandwich factors.
    out=Counter();g=group(f)
    for h,k,sigma,transpose in product(g,g,permutations(range(3)),range(2)):
        w=[int(MUL[MUL[h,j],k]) for j in word]
        if transpose:w=[INV[j] for j in w]
        out[tuple(w[j] for j in sigma)+(w[3],)]+=1
    return out,12*len(g)**2

def projector_dimensions(aux):
    # Trace of an idempotent on Sym^3(A_4) tensor A_2. No characters/irreps.
    dims={}
    for f in ('31','22'):
        total=Fraction(0)
        for h,k,tr in product(group(f),group(f),range(2)):
            action=[int(MUL[MUL[h,INV[j] if tr else j],k]) for j in range(24)]
            traces=[];labels=list(range(24))
            for power in (1,2,3):
                labels=[action[j] for j in labels]
                traces.append(sum(i==j for i,j in enumerate(labels)))
            t1,t2,t3=traces
            ta=sum(aux[i,action[j]] for i,j in enumerate(AUX))
            total+=Fraction(int((t1**3+3*t1*t2+2*t3)*ta),6)
        answer=total/(2*len(group(f))**2)
        assert answer.denominator==1
        dims[f]=int(answer)
    assert dims=={'31':667,'22':1407}
    return dims

def source_audit():
    start=time.monotonic();source=read(ROOT/'source.json')
    aux,conditional,marginals=local_data()
    dims=projector_dimensions(aux)
    output_triples=list(combinations_with_replacement(range(6),3))
    oi={p:i for i,p in enumerate(output_triples)}
    equalities=sparse.load_npz(ROOT/'equality_integer.npz').tocsc()
    offset=0;objective=[];minor_determinants={}
    swap=PI[(3,1,2,0)]
    tp=np.trace(literal(4),axis1=1,axis2=2)
    tq=np.trace(literal(2),axis1=1,axis2=2)
    wp=np.array([np.trace(literal(4)[swap]@a) for a in literal(4)])
    wq=np.array([np.trace(literal(2)[swap]@a) for a in literal(2)])
    for f in FAMILIES:
        B=sparse.load_npz(ROOT/('diagrams_'+f+'.npz')).tocsc()
        words=source['basis_words'][f];bd=source['basis_verification'][f]['diagram_denominator']
        pivotrows=source['basis_verification'][f]['pivot_rows'];columns=[]
        for col,word in enumerate(words):
            terms,den=generator(word,f);assert den==bd
            saved={int(i):int(v) for i,v in zip(B[:,col].indices,B[:,col].data)}
            assert saved=={int(np.ravel_multi_index(w,(24,)*4)):v for w,v in terms.items()}
            # Reduce literal auxiliary permutations and combine physical orbits.
            coeff=Counter()
            for w,n in terms.items():
                base=((min(w[:3])*24+sorted(w[:3])[1])*24+max(w[:3]))*14
                for a in range(14):
                    if aux[a,w[3]]:coeff[base+a]+=n*aux[a,w[3]]
            scale=math.lcm(*(int(v.q) for v in coeff.values()))
            columns.append({i:int(v*scale) for i,v in coeff.items() if v})
            # Independent literal partial traces give each equality coefficient.
            ec=Counter()
            ops=[(conditional[3 if f=='13' else 0],{'31':0,'22':784,'13':1568}[f],1)]
            if f in ('31','22'):ops.append((marginals['aab'],2352,1 if f=='31' else -1))
            if f in ('22','13'):ops.append((marginals['abb'+f],2632,1 if f=='22' else -1))
            collapsed=Counter()
            for w,n in terms.items():collapsed[tuple(sorted(w[:3]))+(w[3],)]+=n
            for (ph,aq),begin,sgn in ops:
                for w,n in collapsed.items():
                    # Physical trace of a permutation has one surviving diagram.
                    choices=[[(j,v) for j,v in enumerate(ph[p]) if v] for p in w[:3]]
                    for choose in product(*choices):
                        j=oi[tuple(sorted(v[0] for v in choose))]
                        v=sgn*n*math.prod(v[1] for v in choose)/den
                        for a,c in enumerate(aq[w[3]]):
                            if c:ec[begin+j*len(aq[w[3]])+a]+=v*c
            observed={int(i):int(v) for i,v in zip(equalities[:,offset].indices,equalities[:,offset].data)}
            expected={i:int(v*source['equality_denominator']) for i,v in ec.items() if v}
            assert all((v*source['equality_denominator']).q==1 for v in ec.values())
            assert observed==expected,(f,col,'equality')
            obj=sum(n*math.prod(int(2*tp[p]-wp[p]) for p in w[:3])*int(wq[w[3]]) for w,n in terms.items())
            obj=Fraction(obj,8*den*128**4) if f=='31' else Fraction(0)
            assert obj==Fraction(source['objective'][offset])
            objective.append(str(obj));offset+=1
        # A nonzero modular minor proves independence over Q. Above idempotent
        # trace proves a matching upper bound, giving complete spanning.
        det=int(nmod_mat([[c.get(p,0) for c in columns] for p in pivotrows],1000033).det())
        assert det
        minor_determinants[f]=det
        print(json.dumps(dict(stage='source',family=f,columns=len(words),seconds=time.monotonic()-start)),flush=True)
    report=dict(source_dimensions=dims,total=2741,basis_minor_prime=1000033,
                nonzero_minor_determinants=minor_determinants,
                source_diagrams_exact=True,equality_columns_exact=2741,
                objective_coefficients_exact=2741,coordinate_scale=128**4,
                all_constraints_rebuilt_from_literal_partial_traces=True)
    (OUT/'source_audit.json').write_text(json.dumps(report,indent=2))
    (OUT/'objective.json').write_text(json.dumps(objective))
    return report

def module_audit():
    cases=read(ROOT/'aligned_modules.json')['representations'];checks=0;count=0
    for case in cases:
        diagrams=literal(case['d'],4,tuple(case['cut']))
        for m in case['modules']:
            J=rational(m['embedding']);G=rational(m['gram'])
            assert J.T*J==G
            # Rank of a literal embedding, not a numeric spectrum.
            assert J.rank()==m['dim']
            for a,r in zip(diagrams,m['basis']):
                assert sp.MutableSparseMatrix(sp.Matrix(a.tolist()))*J==J*rational(r)
                checks+=1
            count+=1
    report=dict(modules=count,literal_intertwiners=checks,positive_metrics=True)
    (OUT/'module_audit.json').write_text(json.dumps(report,indent=2))
    return report

def integer_matrix(rows):
    a=np.array([[Fraction(v) for v in row] for row in rows],dtype=object)
    den=math.lcm(*(v.denominator for v in a.flat))
    return np.array(a*den,dtype=np.int64),den

def maximum(a):
    return max(abs(int(a.min())),abs(int(a.max()))) if a.size else 0

def checked_dot(a,b):
    # Use Python integers to establish the bound, then integer arithmetic only.
    bound=max(sum(abs(int(v)) for v in row) for row in a)*maximum(b)
    assert bound<2**62
    return a@b

CTX=None
def map_init():
    global CTX
    archive=read(ROOT/'aligned_modules.json');lookup={}
    for case in archive['representations']:
        for m in case['modules']:
            G,gd=integer_matrix(m['gram']);r,rd=integer_matrix(np.array(m['basis']).reshape(24*m['dim'],m['dim']))
            h=np.array([checked_dot(G,a) for a in r.reshape(24,m['dim'],m['dim'])])
            lookup[case['d'],tuple(case['cut']),m['label']]=(h,gd*rd,m['dim'])
    blocks=read(ROOT/'positive_blocks.json')['blocks']
    source=read(ROOT/'source.json')
    diagrams={f:sparse.load_npz(ROOT/('diagrams_'+f+'.npz')).T.tocsr() for f in FAMILIES}
    # Independently compose the literal permutation matrices for the joint maps.
    ds=literal(4);tau=PI[(0,1,3,2)];cycle=PI[(2,3,1,0)]
    keys={a.tobytes():i for i,a in enumerate(ds)}
    mappings={'identity':list(range(24)),
              'diag':[keys[(ds[tau]@a@ds[tau].T).tobytes()] for a in ds],
              'cross':[keys[(a@ds[cycle]).tobytes()] for a in ds]}
    CTX=lookup,blocks,source,diagrams,mappings,read(RUN/'exact_frames.json'),read(RUN/'frames.json')['records']

def contract(weight,mods,mapping):
    dims=[m[2] for m in mods]
    value=weight.reshape(dims+dims).transpose(0,4,1,5,2,6,3,7).reshape([d*d for d in dims])
    # Contract directly against G rho, not against rho with a global G weight.
    for axis in sorted(range(4),key=lambda i:dims[i],reverse=True):
        value=np.moveaxis(value,axis,0);shape=value.shape[1:]
        value=checked_dot(mods[axis][0][mapping].reshape(24,-1),value.reshape(value.shape[0],-1)).reshape((24,)+shape)
        value=np.moveaxis(value,0,axis)
    return value.ravel()

def map_one(index):
    lookup,blocks,source,diagrams,mappings,frames,records=CTX
    e=frames[index];rec=records[index];b=blocks[rec['parent']]
    assert e['index']==rec['index']==index and e['parent']==rec['parent']
    if b['type']=='joint':assert b['cut']==[3]
    mods=[lookup[4,tuple(b['cut']),label] for label in b['physical_labels']]+[lookup[2,tuple(b['cut']),b['auxiliary_label']]]
    N,nd=integer_matrix(e['basis']);k=N.shape[1];n=math.prod(m[2] for m in mods)
    assert k==rec['order'] and N.shape[0]==n*(2 if b['type']=='joint' else 1)
    families=FAMILIES if b['type']=='joint' else (b['family'],)
    bd=math.lcm(*(source['basis_verification'][f]['diagram_denominator'] for f in families))
    den=nd**2*math.prod(m[1] for m in mods)*bd
    result=np.zeros((sum(diagrams[f].shape[0] for f in families),k*(k+1)//2),dtype=np.int64)
    for column,(j,i) in enumerate((j,i) for j in range(k) for i in range(j,k)):
        def weight(left,right,factor=1):
            out=np.outer(left[:,i],right[:,j])
            if i!=j:out+=np.outer(left[:,j],right[:,i])
            assert maximum(left)*maximum(right)*factor*2<2**62
            return factor*out
        if b['type']=='source':ops=[(b['family'],weight(N,N),'identity')]
        else:ops=[('31',weight(N[:n],N[:n]),'diag'),('22',weight(N[:n],N[n:],2),'cross'),('13',weight(N[n:],N[n:]),'diag')]
        begin=0
        for f,w,name in ops:
            values=contract(w,mods,mappings[name]);B=diagrams[f]
            bound=int(np.max(abs(B).sum(axis=1)))*maximum(values)
            factor=bd//source['basis_verification'][f]['diagram_denominator']
            assert bound*factor<2**62
            result[begin:begin+B.shape[0],column]=(B@values)*factor;begin+=B.shape[0]
    gcd=math.gcd(den,int(np.gcd.reduce(result.ravel())));result//=gcd;den//=gcd
    saved=np.load(RUN/'rational_maps'/('%04d.npz'%index))
    assert den==int(saved['denominator']) and np.array_equal(result,saved['integer']),('map mismatch',index)
    assert list(families)==saved['families'].tolist()
    np.savez_compressed(OUT/'maps'/('%04d.npz'%index),integer=result,denominator=den,families=families)
    return index

def maps_audit(workers):
    (OUT/'maps').mkdir(exist_ok=True);start=time.monotonic();count=0
    with mp.get_context('spawn').Pool(workers,initializer=map_init) as pool:
        for count,index in enumerate(pool.imap_unordered(map_one,range(1239)),1):
            if count%100==0:print(json.dumps(dict(stage='maps',completed=count,seconds=time.monotonic()-start)),flush=True)
    report=dict(maps_rebuilt=count,all_saved_maps_match=True,arithmetic='bounded integers, no floating point',seconds=time.monotonic()-start)
    (OUT/'maps_audit.json').write_text(json.dumps(report,indent=2));return report

def certificate_audit():
    start=time.monotonic();records=read(RUN/'frames.json')['records']
    target=fmpq_mat([[fmpq(v)] for v in read(OUT/'objective.json')]);total=fmpq_mat(2741,1)
    starts={'31':0,'22':667,'13':2074};rows=0
    for rec in records:
        idx=rec['index'];b=read(CERT/('%04d.json'%idx));k=rec['order']
        assert (b['index'],b['parent'],b['sector'],b['order'])==(idx,rec['parent'],rec['sector'],k)
        Z=fmpq_mat([[fmpq(v) for v in r] for r in b['Z']]);R=fmpq_mat([[fmpq(v) for v in r] for r in b['congruence']])
        assert Z.nrows()==Z.ncols()==R.nrows()==R.ncols()==k and Z==Z.transpose()
        V=R*Z*R.transpose()
        for i in range(k):
            assert V[i,i]-sum(abs(V[i,j]) for j in range(k) if j!=i)>fmpq(9999992,10000000)
            rows+=1
        data=np.load(OUT/'maps'/('%04d.npz'%idx))
        z=fmpq_mat([[Z[i,j]] for j in range(k) for i in range(j,k)])
        v=fmpq_mat(data['integer'].tolist())*z/int(data['denominator'])
        begin=0 if rec['type']=='joint' else starts[rec['family']]
        for i in range(v.nrows()):total[begin+i,0]+=v[i,0]
        if idx%200==0:print(json.dumps(dict(stage='rational certificate',blocks=idx,seconds=time.monotonic()-start)),flush=True)
    eq=read(CERT/'equality.json');E=sparse.load_npz(ROOT/'equality_integer.npz')
    total+=fmpq_mat(E[eq['rows'],:].T.toarray().tolist())*fmpq_mat([[fmpq(v)] for v in eq['coefficients']])
    assert total==target
    report=dict(coefficient_identities=2741,positive_Gram_blocks=1239,positive_congruence_rows=rows,
                rational_margin_lower_bound='9999992/10000000',all_exact=True,
                maps_from_second_implementation=True,seconds=time.monotonic()-start)
    (OUT/'certificate_audit.json').write_text(json.dumps(report,indent=2));return report

def hashes():
    old=read(CERT/'manifest.json');bad=[]
    for name,digest in old.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:bad.append(name)
    assert not bad,bad
    report=dict(existing_manifest_checked=True,files=len(old),all_hashes_match=True)
    (OUT/'manifest_audit.json').write_text(json.dumps(report,indent=2));return report

def reduction_audit():
    # Universal entrywise identities. x/y denote arbitrary multi-indices in V.
    # Monomials are represented as (vector, multi-index, conjugated) factors.
    x=['x'+str(i) for i in range(4)];y=['y'+str(i) for i in range(4)]
    row=x[:3]+[y[3]];col=y[:3]+[x[3]]  # Last-replica transpose.
    tau=(0,1,3,2);cycle=(2,3,1,0)
    def permute(v,p):return [v[p.index(i)] for i in range(4)]
    def density(v,r,c):return Counter([(v[i],r[i],False) for i in range(4)]+[(v[i],c[i],True) for i in range(4)])
    z=[('a',x[0],False),('a',x[1],False),('b',x[2],False),('a',x[3],True)]
    wbar=[('b',y[0],True),('b',y[1],True),('a',y[2],True),('b',y[3],False)]
    assert density('aabb',row,permute(col,cycle))==Counter(z+wbar)
    for v in ('aaab','bbba'):
        first='a' if v=='aaab' else 'b';second='b' if first=='a' else 'a'
        expected=[(first,x[0],False),(first,x[1],False),(second,x[2],False),(first,x[3],True),
                  (first,y[0],True),(first,y[1],True),(second,y[2],True),(first,y[3],False)]
        assert density(v,permute(row,tau),permute(col,tau))==Counter(expected)
    # This part is a diagnostic on general complex inputs, not a proof tolerance.
    rng=np.random.default_rng(202610140731);worst=0.;unequal=0.;rank_one=0
    for test in range(24):
        A=np.linalg.qr(rng.normal(size=(64,2))+1j*rng.normal(size=(64,2)))[0]
        B=rng.normal(size=(64,2))+1j*rng.normal(size=(64,2))
        if test%3==0:B[:,1]=0;rank_one+=1
        if test%3==1:B[:,1]*=.13
        B*=np.sqrt(2/np.vdot(B,B).real)
        C=A@B.conj().T/2
        assert np.linalg.norm(C-A@(2*C.conj().T@A).conj().T/2)<1e-12
        a=A.reshape(4,4,4,2)/np.sqrt(2);b=B.reshape(4,4,4,2)/np.sqrt(2)
        psi=np.multiply.outer(a,b)
        q=0.;swap_value=0.;werner_value=0.
        for mask in range(8):
            sites=[i for i in range(3) if mask>>i&1]
            D=C.reshape((4,)*6);remaining=3
            for i in reversed(sites):D=np.trace(D,axis1=i,axis2=remaining+i);remaining-=1
            term=np.vdot(D,D).real
            order=list(range(8))
            for i in sites+[3]:order[i],order[i+4]=order[i+4],order[i]
            observed=np.vdot(psi,psi.transpose(order))
            worst=max(worst,abs(observed-term))
            q+=(-.5)**len(sites)*term;swap_value+=(-.5)**len(sites)*observed
            # Literal Omega_i action on vec(C), with Omega=sum_j |jj>.
            v=C.reshape((4,)*6)
            for i in sites:
                tr=np.trace(v,axis1=i,axis2=i+3)
                w=np.zeros_like(v)
                for j in range(4):
                    idx=[slice(None)]*6;idx[i]=idx[i+3]=j;w[tuple(idx)]=tr
                v=w
            werner_value+=(-.5)**len(sites)*np.vdot(C.ravel(),v.ravel())
        worst=max(worst,abs(q-swap_value),abs(q-werner_value))
        unequal=max(unequal,float(np.linalg.norm(B.conj().T@B-np.eye(2))))
    assert worst<1e-12
    report=dict(universal_joint_entrywise_identities=3,
                general_complex_d4_diagnostic_inputs=24,rank_one_inputs=rank_one,
                swap_subset_diagnostics=192,maximum_diagnostic_error=float(worst),
                maximum_B_isometry_defect=unequal,numerical_checks_are_not_proof=True,
                analytic_reduction='See manuscript, matched four-replica section')
    (OUT/'reduction_audit.json').write_text(json.dumps(report,indent=2));return report

def equality_rank_audit():
    start=time.monotonic();E=sparse.load_npz(ROOT/'equality_integer.npz')
    rows=read(CERT/'equality.json')['rows']
    rank=fmpq_mat(E.toarray().tolist()).rank()
    selected=fmpq_mat(E[rows,:].toarray().tolist()).rank()
    assert rank==selected==272
    report=dict(all_row_rank=rank,selected_row_rank=selected,arithmetic='rational',seconds=time.monotonic()-start)
    (OUT/'equality_rank_audit.json').write_text(json.dumps(report,indent=2));return report

def main():
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['foundation','maps','certificate','reduction','rank']);parser.add_argument('--workers',type=int,default=3);args=parser.parse_args()
    OUT.mkdir(exist_ok=True)
    if args.stage=='foundation':
        print(json.dumps(hashes()),flush=True)
        print(json.dumps(module_audit()),flush=True)
        print(json.dumps(source_audit()),flush=True)
    elif args.stage=='maps':print(json.dumps(maps_audit(args.workers)),flush=True)
    elif args.stage=='reduction':print(json.dumps(reduction_audit()),flush=True)
    elif args.stage=='rank':print(json.dumps(equality_rank_audit()),flush=True)
    else:print(json.dumps(certificate_audit()),flush=True)

if __name__=='__main__':main()
