classdef InitialCondition
    %INITIALCONDITION Summary of this class goes here
    %   Detailed explanation goes here

    properties(Constant)
        EXPONENTIAL = 0;
        GAUSSIAN = 1;
        PLANE = 2;
        SECH = 3;
    end

    properties
        n;
        SELECTION = "";
    end

    methods
        function obj = InitialCondition(selection, n_, A0_, x_, y_)
            switch(selection)
                case obj.EXPONENTIAL
                    obj.n = obj.initial_exp(n_, A0_, x_, y_);
                    obj.SELECTION = "EXPONENTIAL";
                case obj.GAUSSIAN
                    obj.n = obj.initial_gaussian_2d(n_, A0_);
                    obj.SELECTION = "GAUSSIAN";
                case obj.PLANE
                    obj.n = obj.initial_plane_wave(n_, A0_);
                    obj.SELECTION = "PLANE";
                case obj.SECH
                    obj.n = obj.initial_sech(n_, A0_, 0);
                    obj.SELECTION = "SECH";
                otherwise
                    obj.n = null;
            end
        end

        function n = initial_exp(obj, n, A0, x, y)
            dims = size(n);
            
            X = x .* ones(dims(1), dims(2));
            Y = y' .* ones(dims(1), dims(2));
            
            n = A0 * exp(-2*((X - 3.75).^2 + (Y - 3.75).^2));
        end

        function n = initial_gaussian_2d(obj, n, A0)
            dims = size(n);
            a = A0;
            b1 = round(dims(2)/2);
            b2 = round(dims(1)/2);
            c = sqrt(dims(1));
            
            for k=1:dims(1)
                for j=1:dims(2)
                    n(k,j) = n(k,j) + a * exp(-((j-b1).^2 + (k-b2).^2) / (2*c^2));
                end
            end
        end

        function n = initial_plane_wave(obj, n, A0)
            dims = size(n);
            
            t = 0:dims(2)-1;
            a = A0;
            b1 = round(dims(2)/3);
            c = 100/dims(2);
            
            for i=1:dims(1)
                n(i,:) = n(i,:) + a * exp(-(t-b1).^2 / (2*c^2));
            end
        end

        % Remember to set the constants in this function according to Wei1995 when
        % doing actual final tests.
        % Set is_velocity to 1 if you want to get the velocity function instead
        function n = initial_sech(obj, n, A0, is_velocity)
            dims = size(n);
            
            x = 0:dims(2)-1;
            a1 = A0/2;
            a2 = A0/2;
            B = 0.2;
            b1 = dims(2)/4;
            c = 1;
            
            smaller_dim = min(dims(2),dims(1));
            
            if is_velocity == 0
                for i=1:smaller_dim
                    n(i,:) = n(i,:) + a1 * sech(B*(c*x-b1)).^2 + a2 * sech(B*(c*x-b1)).^4;
                end
            elseif(is_velocity == 1)
                for i=1:smaller_dim
                    n(i,:) = n(i,:) + a1 * sech(B*(c*x-b1)).^2;
                end
            end
        end
    end

    methods(Static)
        function profileEnum = getProfileEnum(profileName)
            switch(profileName)
                case "EXPONENTIAL"
                    profileEnum = InitialCondition.EXPONENTIAL;
                case "GAUSSIAN"
                    profileEnum = InitialCondition.GAUSSIAN;
                case "PLANE"
                    profileEnum = InitialCondition.PLANE;
                case "SECH"
                    profileEnum = InitialCondition.SECH;
                otherwise
                    error('You have provided an invalid InitialCondition name.');
            end
        end
    end
end

