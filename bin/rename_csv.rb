#!/usr/bin/env ruby

#name=ARGV.shift
sample_list=ARGV.shift
samples=File.readlines(sample_list)

cmd="mkdir exome_depth_results"
system(cmd)
Dir.glob('*.csv').each do |name|

number=name.split('_').last.split('.').first.to_i
cmd = "cp #{name} exome_depth_results/#{samples[number].gsub('exome_sorted.bam','cnv.csv')}"
system(cmd)
end
