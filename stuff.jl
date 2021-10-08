using SDDP

abstract type Test end

print(Test)
print(typeof(Test))

mutable struct Subtest <: Test end

print(Subtest)
print(typeof(Subtest))
print(typeof(Subtest) == Subtest)
print(Subtest == Subtest)
lim = SDDP.IterationLimit
print(lim)
print(typeof(lim))
print(lim == SDDP.IterationLimit)
