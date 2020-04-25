#-------------------------------------------------------------------------------
## BRENT'S METHOD ##
#-------------------------------------------------------------------------------
function Brent(f::Function, a::Real, b::Real; tol = 1e-14)
    # Implementation of Brent's method to find a root of a function (as on wikipedia)
    fa = f(a); fb = f(b)
    if fa * fb > 0
        error("f[a] and f[b] should not have different signs!")
    end

    c = a;
    fc=fa;   # at the beginning: c = a
    c = a;
    d = b - a;
    e = d

    iter = 0
    maxiter = 10000

    while iter < maxiter
        iter += 1

        if fb * fc > 0
            c = a; fc = fa; d = b - a; e = d
        end

        if abs(fc) < abs(fb)
            a = b; b = c; c = a
            fa = fb; fb = fc; fc = fa
        end
        tol = 2.0 * eps() * abs(b) + tol; m = (c - b) / 2.0; #Toleranz

        if abs(m)>tol && abs(fb)>0 #Verfahren muss noch durchgeführt werden
            if abs(e)<tol || abs(fa)<=abs(fb)
                d=m; e=m
            else
                s=fb/fa
                if a==c
                    p=2*m*s; q=1-s
                else
                    q=fa/fc; r=fb/fc
                    p=s*(2*m*q*(q-r)-(b-a)*(r-1))
                    q=(q-1)*(r-1)*(s-1)
                end
                if p>0
                    q=-q
                else
                    p=-p
                end
                s=e; e=d
                if  2*p<3*m*q-abs(tol*q)  && (p<abs(s*q/2))
                    d=p/q
                else
                    d=m; e=m
                end
            end
            a=b; fa=fb
            if abs(d)>tol
                b=b+d
            else
                if m>0
                    b=b+tol
                else
                    b=b-tol
                end
            end
        else
            break
        end
        fb = f(b)
    end
    return b, iter
end
