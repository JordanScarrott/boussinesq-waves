% Creates a vector field of size dims = (rows cols)
classdef Vector_Field
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        u
        v
    end
    
    
    methods
    % Constructor
        function obj = Vector_Field(dims)
            obj.u = zeros(dims(1), dims(2));
            obj.v = zeros(dims(1), dims(2));
        end
    end
    
    % Static Methods
    methods(Static)
        % Static add - Adds two vector fields element wise
        function obj = add(vecf1, vecf2)
            dims = size(vecf1);
            obj = Vector_Field([dims(1), dims(2)]);
            obj.u = vecf1.u + vecf2.u;
            obj.v = vecf1.v + vecf2.v;
        end
        
    end
    
end

