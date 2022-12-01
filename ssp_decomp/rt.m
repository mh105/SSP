function rtB= rt(B)
assert(all(size(B)==[2,2]))
rtB= B(2,1)-B(1,2);
end