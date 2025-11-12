classdef DataStorage < handle  %  <===== THIS LINE MAKES IT A HANDLE CLASS
    properties 
        % We can store any amount of data in this empty structure
        dataArea = struct(); 
    end
    
    methods 
        function obj = DataStorage(obj)
            return;
        end
    end
end