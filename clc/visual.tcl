set molnum 0

mol load pdb $argv
mol delrep 0 $molnum

mol representation NewCartoon
mol color Structure
mol addrep $molnum
