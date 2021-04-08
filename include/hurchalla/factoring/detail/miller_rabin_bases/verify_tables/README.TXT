I generated all the miller-rabin hash tables myself (via brute-force searches, the details of which are beyond the scope of discussion here).
I verified their correctness using the verification scripts/cpp files in this folder.  If you wish to verify the hash tables and bases yourself, you can view and run the scripts/cpp files in this folder.  You will need a file that lists the strong base2 psuedoprimes for all numbers < 2^64:

You can get a two-part 7zip archive that contains the desired file, "strong_psps_to_2_64.txt", from
https://raw.githubusercontent.com/hurchalla/Extra-Files/main/psps.7z.001
https://raw.githubusercontent.com/hurchalla/Extra-Files/main/psps.7z.002

( I programmatically created the "strong_psps_to_2_64.txt" file using Jan Feitsma's Fermat base2 pseudoprimes for all numbers < 2^64, available at http://www.cecm.sfu.ca/Pseudoprimes/index-2-to-64.html  and see also http://www.janfeitsma.nl/math/psp2/index )
