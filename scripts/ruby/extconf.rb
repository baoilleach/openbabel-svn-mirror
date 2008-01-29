# Creates a makefile to build the Open Babel-Ruby extension.

# Compensate for the fact that OpenBabel on OS X will not compile
# as a universal binary by default.
if system('uname > /dev/null') && `uname -smr` =~ /^Darwin 9.* i386$/
  ENV['ARCHFLAGS'] = '-arch i386'
end

require 'mkmf'

dir_config('openbabel')

# Find a trivial header in order to add the proper include path
# to the build flags.
find_header('inchi_api.h', '../../include')

# Prevent Ruby 1.8.x from trying to compile and link the extension
# using gcc.
if RUBY_VERSION < "1.9"
  cxx = ''
  begin
    File.open('../Makefile', 'r').each_line do |line|
      if line =~ /CXX = /
        cxx = Regexp.last_match.post_match.chomp
      end
    end
  rescue Errno::ENOENT
    puts 'Please configure OpenBabel before compiling the Ruby extension'
  end
  
  cpp_command(cxx)
  CONFIG['LDSHARED'] = "#{cxx} #{ENV['ARCHFLAGS']} -pipe -bundle"
end

if have_library('openbabel')
  create_makefile('openbabel')
else
  puts "Install Open Babel first. If you've already compiled and installed Open Babel, you may need to run ldconfig."
end

