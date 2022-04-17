classdef FloorProfile
    %FLOORPROFILE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
%         Here we speak in terms of the floor height above zero instead of
%         the depth of the water but what is returned is the water depth
        FLAT = [[0, 0]; [0, 0]; [0, 0]; [0, 0]];
        SINGLE_BAR = [[0, 0]; [0.25, 0]; [0.35, 1/3]; [0.4, 1/3]; [0.45, 1/5]; [0.5, 1/5]; [0.75, 1]; [1, 1]];
    end

    properties
        floorProfile;
        y_data;
    end
    
    methods
        function obj = FloorProfile(preset, x_, h0_)
            obj.floorProfile = Profile_description(preset, x_);
            obj.y_data = obj.floorProfile.y .* h0_;
        end
        
    end
end

