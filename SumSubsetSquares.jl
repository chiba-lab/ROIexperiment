#* Signal processing library 
using DSP
import Base.Threads.@spawn

#* Generating function for distribution of sums:  f(s) = ΣG(x,y) = Σg_{n,k}*x^n*y^k = ∏_n_i(1+x*y^k_i) 
#* where and g_n,k is the number of subsets of size n with k sum  

#* creates vector of factors of generator function P[k] = (1+x*y^k) for sequence of squares from 1:k_max²
function generating_factors(k_max)
    #* n+1 x k+1 matrix encoding: row = size of subset, column = degree of value of sum
    #* first row and column represent n,k=0 
    function f̨ₖ(k)
        f=zeros(Int8,k^2+1,2)
        f[1,1]=1
        f[k^2+1,2]=1
        return f
    end
    return [f̨ₖ(k) for k in 1:k_max]
end

#* Clips values for to prevent overflow
function clipvals(x)
    x=abs(round(x))
    if x==0
        return 0
    elseif x==1
        return 1
    else
        return 2
    end
end

#* Recursively reduce factors by convolution  
function dp_conv_reduce(G)
    s=size(G,1)
    if s==1
        return clipvals.(G[1])
    end
    if s==2
        return clipvals.(conv(G[1], G[2]))
    end
    g1 = @spawn dp_conv_reduce(G[1:Int(floor(s/2))])
    g2 = @spawn  dp_conv_reduce(G[Int(floor(s/2))+1:end])
    return clipvals.(conv(fetch(g1) , fetch(g2)))
end

#* computes coefficients g_n_k for subsets of size n
function compute_generating_coefs(G, n; dp=true)
   #* computes product using convolutions. Convolutions use FFTs to make it go way faster!
    if dp
        coefs= dp_conv_reduce(G)
    else 
        coefs = reduce((x,y)->clipvals.(conv(x,y)), G)
    end
    sub_sums =  collect(0:size(coefs,1)-1)
    return coefs[:, n+1], sub_sums 
end


#* sum of the unique subset sums 
function sum_unique_subset_sums(k_max, n)
    G = generating_factors(k_max)
    coefs, sub_sums = compute_generating_coefs(G, n)
    return sum(sub_sums[coefs.==1])
end

sum_unique_subset_sums(100,50)   