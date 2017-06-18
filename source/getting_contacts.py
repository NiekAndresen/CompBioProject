import contact_tools_fromMahmoud as ct

nat_cst = ct.load_native_contacts("/home/niek/HSA_data/1ao6/1ao6A.pdb")
ct.write_rosetta_cst_restraints(nat_cst, "/home/niek/HSA_data/1ao6/1ao6A.cst", "/home/niek/HSA_data/1ao6/1ao6A.fasta")
