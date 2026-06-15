import streamlit as st
from gene_sequence_analysis import (
    gc_content,
    reverse_complement,
    translate_dna,
    find_orfs,
)

st.title("Gene Sequence Analysis Tool")

seq = st.text_area("Enter DNA Sequence")

if st.button("Analyze"):
    gc = gc_content(seq)

    st.write("GC Content:", f"{gc:.2f}%")
    st.write("Reverse Complement:")
    st.code(reverse_complement(seq))

    st.write("Protein Translation:")
    st.code(translate_dna(seq))

    orfs = find_orfs(seq)

    st.write(f"ORFs Found: {len(orfs)}")

    for i, orf in enumerate(orfs[:5]):
        st.write(
            f"ORF {i+1}: {orf['start']} - {orf['end']} "
            f"({orf['length']} bp)"
        )
