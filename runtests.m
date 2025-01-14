addpath('tests/')
tests = dir(fullfile('tests/','*.m'));
pass_count = 0;
test_count = 0;
for k = 1:length(tests)
  test = tests(k).name;
  run(test);
  passed = ans;
  pass_count = pass_count + passed;
  test_count = test_count + 1;
end

if (pass_count < test_count)
    warning("\n\n%i out of %i tests passed\n", pass_count, test_count);
else
    fprintf("\n\nAll %i tests passed\n", test_count);
end
