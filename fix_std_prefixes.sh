#!/bin/bash
# Script to fix missing std:: prefixes in MOLSET source files
# This fixes common C++ standard library types/functions that need std:: namespace

MOLSET_DIR="/srv/disk1/MYPROG/HAPACK/MOLSET/sources"

echo "Fixing std:: prefixes in MOLSET source files..."

# Process all .cpp and .h files
for file in "$MOLSET_DIR"/*.cpp "$MOLSET_DIR"/*.h; do
    if [ -f "$file" ]; then
        echo "Processing: $(basename "$file")"

        # Create backup
        cp "$file" "$file.bak"

        # Fix cout, cerr, cin (but not if already has std::)
        sed -i 's/\([^:a-zA-Z_.]\)cout\s*<</\1std::cout <</g' "$file"
        sed -i 's/\([^:a-zA-Z_.]\)cerr\s*<</\1std::cerr <</g' "$file"
        sed -i 's/\([^:a-zA-Z_.]\)cin\s*>>/\1std::cin >>/g' "$file"

        # Fix endl (but not if already has std::)
        sed -i 's/<<\s*endl\s*;/<< std::endl;/g' "$file"
        sed -i 's/<<\s*endl\s*$/<< std::endl/g' "$file"

        # Fix stream types at start of declarations (with tab or spaces)
        sed -i 's/^\(\s*\)ifstream\s\+/\1std::ifstream /g' "$file"
        sed -i 's/^\(\s*\)ofstream\s\+/\1std::ofstream /g' "$file"
        sed -i 's/^\(\s*\)fstream\s\+/\1std::fstream /g' "$file"
        sed -i 's/^\(\s*\)istrstream\s\+/\1std::istrstream /g' "$file"
        sed -i 's/^\(\s*\)ostrstream\s\+/\1std::ostrstream /g' "$file"
        sed -i 's/^\(\s*\)stringstream\s\+/\1std::stringstream /g' "$file"
        sed -i 's/^\(\s*\)istringstream\s\+/\1std::istringstream /g' "$file"
        sed -i 's/^\(\s*\)ostringstream\s\+/\1std::ostringstream /g' "$file"

        # Fix stream types in function parameters (type& or type *)
        sed -i 's/(\s*ifstream\s*&/(std::ifstream\&/g' "$file"
        sed -i 's/(\s*ofstream\s*&/(std::ofstream\&/g' "$file"
        sed -i 's/(\s*istream\s*&/(std::istream\&/g' "$file"
        sed -i 's/(\s*ostream\s*&/(std::ostream\&/g' "$file"
        sed -i 's/,\s*ifstream\s*&/, std::ifstream\&/g' "$file"
        sed -i 's/,\s*ofstream\s*&/, std::ofstream\&/g' "$file"
        sed -i 's/,\s*istream\s*&/, std::istream\&/g' "$file"
        sed -i 's/,\s*ostream\s*&/, std::ostream\&/g' "$file"

        # Fix standalone istream/ostream in function signatures
        sed -i 's/^\([^/]*\)ostream\s*&\s*\([a-zA-Z_]\)/\1std::ostream\& \2/g' "$file"
        sed -i 's/^\([^/]*\)istream\s*&\s*\([a-zA-Z_]\)/\1std::istream\& \2/g' "$file"

        # Fix vector declarations (careful not to match HaVec or other *Vec*)
        sed -i 's/^\(\s*\)vector\s*</\1std::vector</g' "$file"
        sed -i 's/\(\s\)vector\s*</\1std::vector</g' "$file"
        sed -i 's/^\(\s*\)vector\s*::/\1std::vector::/g' "$file"

        # Fix list declarations
        sed -i 's/^\(\s*\)list\s*</\1std::list</g' "$file"
        sed -i 's/\(\s\)list\s*</\1std::list</g' "$file"

        # Fix set declarations (careful - "set" is common word)
        sed -i 's/^\(\s*\)set\s*</\1std::set</g' "$file"
        sed -i 's/\(\s\)set\s*</\1std::set</g' "$file"

        # Fix map declarations
        sed -i 's/^\(\s*\)map\s*</\1std::map</g' "$file"
        sed -i 's/\(\s\)map\s*</\1std::map</g' "$file"

        # Fix pair
        sed -i 's/^\(\s*\)pair\s*</\1std::pair</g' "$file"
        sed -i 's/\(\s\)pair\s*</\1std::pair</g' "$file"

        # Fix less<> comparator
        sed -i 's/,\s*less\s*</, std::less</g' "$file"
        sed -i 's/<\s*less\s*</< std::less</g' "$file"

        # Fix string (careful - very common, only fix obvious cases)
        sed -i 's/^\(\s*\)string\s\+\([a-zA-Z_][a-zA-Z0-9_]*\)\s*[;=]/\1std::string \2;/g' "$file"
        sed -i 's/^\(\s*\)string\s*&/\1std::string\&/g' "$file"
        sed -i 's/(\s*string\s*&/(std::string\&/g' "$file"
        sed -i 's/,\s*string\s*&/, std::string\&/g' "$file"

        # Fix make_pair
        sed -i 's/\([^:a-zA-Z_.]\)make_pair\s*(/\1std::make_pair(/g' "$file"

        # Fix iterator
        sed -i 's/>::iterator/>::iterator/g' "$file"  # Already correct, just ensure consistency

        # Fix getline (but not member function .getline())
        sed -i 's/\([^:a-zA-Z_.]\)getline\s*(/\1std::getline(/g' "$file"

        # Fix sort, min, max (common algorithms - skip copy/find/swap as they may be member functions)
        sed -i 's/\([^:a-zA-Z_.]\)sort\s*(/\1std::sort(/g' "$file"
        sed -i 's/\([^:a-zA-Z_.]\)min\s*(/\1std::min(/g' "$file"
        sed -i 's/\([^:a-zA-Z_.]\)max\s*(/\1std::max(/g' "$file"

        # Clean up any double std::std:: that might have been introduced
        sed -i 's/std::std::/std::/g' "$file"

    fi
done

echo "Done! Backup files created with .bak extension"
echo "To remove backups: rm $MOLSET_DIR/*.bak"
