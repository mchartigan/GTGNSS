classdef EmptyMeasurement < Measurement
    %EMPTYMEASUREMENT Provides a default instantiation of abstract class
    %Measurement.

    properties
        dim = 1
    end

    methods
        function y = computemeas(~,~,~)
            y = 0;
        end

        function H = measpartials(~,~,~)
            H = 0;
        end
    end
end