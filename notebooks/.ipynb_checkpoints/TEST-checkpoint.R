# max.R
# Fetch command line arguments
# Convert to numerics
nums = c(5,6,7,8,9,2)
# cat will write the result to the stdout stream
cat(max(nums))
print('hi')

write.table(nums,"test_nums.txt",sep="\t",row.names=FALSE)