% AnotherExample.m
classdef AnotherExample < Example % AnotherExample inherits from Example
  properties
    property3 = 5
  end
  methods
    function obj = AnotherExample(input1)
      obj.property3 = input1;
    end
    function obj = method2(obj,input1)
      % Do something
      obj.property1 = input1;
    end
    
  end
  methods(Static)
    function main
      Example.main
      fprintf('Constructing AnotherExample\n')
      AE = AnotherExample(.001)
      fprintf('Using method1\n')
      AE = AE.method1(1.1)
      fprintf('Using method2\n')
      AE = AE.method2(2.2)
    end
  end
end

% Put this "classdef" a file called AnotherExample.m 
% Make sure to have Example.m in the same directory.
% Then run
% >> Example.main
