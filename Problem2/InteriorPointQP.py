import numpy as np
from scipy.linalg import ldl, solve_triangular
from numpy.linalg import norm

"""
Primal-Dual Predictor-Corrector Interior Point Quadratic Programming algorithm

min 0.5 x' H x + g' x
subject to A' x  = b
           C' x >= d 

x: variable we want to find
y: equality constraint lagrange multiplier
z: inequality constraint lagrange multiplier
s: slack variable

The algorithm below requires an initial guess
(x0,y0,z0,s0) where z0 > 0 and s0 > 0.
x0 should be an interior point, it must not be on the boundary.


"""


def InteriorPointQP(H,g,A,b,C,d,x0,y0,z0,s0,MaxIter = 100, tol = 10**(-6)):
    
    x,y,z,s = x0,y0,z0,s0
    
    def simple_stopping(r_L,r_A,r_C,mu,tol):
        r = np.block([r_L,r_A,r_C,mu])
        return norm(r,np.inf) < tol
    
    # Calculate residuals
    r_L  = H @ x + g - A @ y - C @ z
    r_A  = b - A.T @ x
    r_C  = s + d - C.T @ x
    r_sz = s*z
    
    n,mc = C.shape
    mu = (z.T @ s)/mc
    
    converged = False
    k = 0
    
    X = np.zeros([MaxIter+1,len(x)])
    X[0,:] = x
    
    while converged == False and k < MaxIter:
        
        # Compute H_bar and setup KKT system
        H_bar = H + C @ np.diag(z/s) @ C.T
        m = A.shape[1]
        KKT = np.block([[H_bar, -A],[-A.T, np.zeros([m,m])]])
        
        # we find ldl factorization of the KKT system
        L, D, perm = ldl(KKT)
        
        # Compute affine direction
        r_L_bar = r_L - C @ np.diag(z/s) @ (r_C - r_sz/z)
        rhs = (-1)*np.block([r_L_bar, r_A])
        rhs2 = solve_triangular(L[perm,:], rhs[perm],lower=True)
        res = solve_triangular(D @ L[perm,:].T, rhs2)[perm]
        dx_aff = res[:len(x)]
        dy_aff = res[len(x):]
        
        dz_aff = (-1)*np.diag(z/s) @ C.T @ dx_aff + np.diag(z/s) @ (r_C - r_sz/z)
        ds_aff = - r_sz/z - (s * dz_aff)/z
        
        # Find larges affine step alpha
        alphas = np.linspace(0,1,50)
        alphas = alphas[:, np.newaxis]
        candidates_aff = np.block([z,s]) + alphas * np.block([dz_aff, ds_aff])
        largest_index_aff = np.argmin(np.all(candidates_aff >= 0, axis=1))
        alpha_aff = alphas[largest_index_aff-1]
        
        # Duality gap and centering parameter
        mu_aff = ((z + alpha_aff * dz_aff).T @ (s + alpha_aff * ds_aff)) / mc
        sigma = (mu_aff/mu)**3
        
        # Affine-Centering-Correction Direction
        r_sz_bar = r_sz + ds_aff*dz_aff - sigma*mu * np.ones(len(r_sz))
        r_L_bar = r_L - C @ np.diag(z/s) @ (r_C - r_sz_bar/z)
        
        rhs = (-1)*np.block([r_L_bar, r_A])
        rhs2 = solve_triangular(L[perm,:], rhs[perm],lower=True)
        res = solve_triangular(D @ L[perm,:].T, rhs2)[perm]
        dx = res[:len(x)]
        dy = res[len(x):]
        
        
        dz = (-1)*np.diag(z/s) @ C.T @ dx + np.diag(z/s) @ (r_C - r_sz_bar/z)
        ds = -r_sz_bar/z - s * dz/z
        
        candidates = np.block([z,s]) + alphas * np.block([dz, ds])
        largest_index = np.argmin(np.all(candidates >= 0, axis=1))
        alpha = alphas[largest_index-1]
    
        
        # Update iterate
        nu = 0.9
        alpha_bar = nu*alpha 
        
        x += alpha_bar * dx
        y += alpha_bar * dy
        z += alpha_bar * dz
        s += alpha_bar * ds
        
        # Calculate residuals
        r_L  = H @ x + g - A @ y - C @ z
        r_A  = b - A.T @ x
        r_C  = s + d - C.T @ x
        r_sz = s*z
        
        mu = (z.T @ s)/mc
        
        converged = simple_stopping(r_L,r_A,r_C,mu,tol)
            
        k += 1
        X[k,:] = x
        
    X = X[:(k+1),:]
    results = dict()
    results["xmin"] = x
    results["lagrange_eq"] = y
    results["lagrange_ineq"] = z
    results["converged"] = converged
    results["iter"] = k
    results["x_array"] = X
    
    return results


