% Example.m
classdef Example
  properties % These can be made private/public/static etc.
    property1     % Default value will be the empty array
    property2 = 3 % Default values can be set
  end
  methods % These can be made private/public/static etc.
    function obj = Example(input1,input2,varargin)
      % The constructor
      if nargin>0 % For inheritance must allow constructing with no inputs
      	obj.property1 = input1*input2;
      	if nargin>2
      	  obj.property2 = varargin{1};
      	end 
      end
    end
    function obj = method1(obj,input1)
      obj.property2 = input1;
    end 
  end
  methods(Static)
    function main
      fprintf('Constructing Example with no inputs\n');
      E = Example
      fprintf('Constructing Example with inputs 10,2,13\n');
      E = Example(10,2,13)
    end
  end
end

% To run this example, put this "classdef" a file named Example.m and run
% >> Example.main
