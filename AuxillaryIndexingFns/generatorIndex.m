function index = generatorIndex(dim, i)

        index = 0;
        if i>1
            for j=1:(i-1)
                index = index + dim;
            end
        end
        
end