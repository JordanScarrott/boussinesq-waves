classdef FloorProfile
    %FLOORPROFILE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        FLAT = 0;
        SINGLE_BAR = 1;

%         Here we speak in terms of the floor height above zero instead of
%         the depth of the water but what is returned is the water depth
        FLAT_DATA = [[0, 0]; [0, 0]; [0, 0]; [0, 0]];
        SINGLE_BAR_DATA = [[0, 0]; [0.25, 0]; [0.35, 1/3]; [0.4, 1/3]; [0.45, 1/5]; [0.5, 1/5]; [0.75, 1]; [1, 1]];
    %             obj.h = obj.h0-obj.h0/(obj.real_x*12.7/13) * ((ones(obj.yn, obj.xn) .* obj.x ./ obj.scale) - (obj.real_x*0.3/13));
    end

    properties
        SELECTION = "";
        profileDescription;
        y_data;
    end
    
    methods
        function obj = FloorProfile(preset, x_, h0_)
            switch(preset)
                case obj.FLAT
                    obj.profileDescription = Profile_description(obj.FLAT_DATA, x_);
                    obj.SELECTION = "FLAT";
                case obj.SINGLE_BAR
                    obj.profileDescription = Profile_description(obj.SINGLE_BAR_DATA, x_);
                    obj.SELECTION = "SINGLE_BAR";
            end
            obj.y_data = obj.profileDescription.y .* h0_;
        end
        
    end
end

