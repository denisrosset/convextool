classdef CVXTestCase < matlab.unittest.TestCase
       methods(TestMethodTeardown)
        function clear(testCase)
            cvx_clear;
        end            
    end
    methods(TestMethodSetup)
        function setSolver(testCase)
            cvx_clear;
            cvx_solver(ConfigTestOption('solver'));
        end
    end 
end
