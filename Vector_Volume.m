% Returns a list of Vector_Fields of depth, iterations.
% Dimensions is an array of the dims of the Vector_Fields (rows, cols)
function S = Vector_Volume(dimensions, iterations)
    
    for i=1:iterations
        S(i) = Vector_Field([dimensions(1), dimensions(2)]);
    end
    
    
end

