var myHeading = document.getElementById("demo");

var text = 'TAATGCCATGGGATGTT';
var number = 3;

function cyclic_spectrum(peptide){
  //   Generate the theoretical spectrum of a cyclic peptide.
  //   Input: An amino acid string Peptide.
  //   Output: Cyclospectrum(Peptide).
  var amino_acid_mass = { 'G' : 57,  'A' : 71,  'S' : 87,  'P' : 97, 'V' : 99,  'T' : 101, 'C' : 103, 'I' : 113, 'L' : 113, 'N' : 114,
                      'D' : 115, 'K' : 128, 'Q' : 128, 'E' : 129,'M' : 131, 'H' : 137, 'F' : 147, 'R' : 156, 'Y' : 163, 'W' : 186};
  var prefix_mass = [0];
  var peptide_length = peptide.length;
  
  for (i = 1; i <= peptide_length; i++){
    prefix_mass.push(prefix_mass[i-1] + amino_acid_mass[peptide[i-1]]);
  }
  
  var total_peptide_mass = prefix_mass[peptide_length];
  
  for (i = 1; i <= peptide_length; i++){
    for (j = i + 1; j <= peptide.length; j++){
      prefix_mass.push(prefix_mass[j] - prefix_mass[i]);
      if (i > 0 && j < peptide_length){
        prefix_mass.push(total_peptide_mass - (prefix_mass[j] - prefix_mass[i]));
      }
    }
  }
  
  prefix_mass.sort(function(a, b) { return a - b; });
  
  return prefix_mass;
}

function rna_to_aa(rna){
  /*Translate an RNA string into an amino acid string.
    Input: An RNA string Pattern.
    Output: The translation of Pattern into an amino acid string Peptide.*/
  
  var rna_codon_table = {'UUU':'F', 'UUC':'F', 'UUA':'L', 'UUG':'L', 'UCU':'S',
                         'UCC':'S', 'UCA':'S', 'UCG':'S', 'UAU':'Y', 'UAC':'Y', 
                         'UAA':'*', 'UAG':'*', 'UGU':'C', 'UGC':'C', 'UGA':'*', 
                         'UGG':'W', 'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L', 
                         'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAU':'H', 
                         'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGU':'R', 'CGC':'R', 
                         'CGA':'R', 'CGG':'R', 'AUU':'I', 'AUC':'I', 'AUA':'I', 
                         'AUG':'M', 'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 
                         'AAU':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGU':'S', 
                         'AGC':'S', 'AGA':'R', 'AGG':'R', 'GUU':'V', 'GUC':'V', 
                         'GUA':'V', 'GUG':'V', 'GCU':'A', 'GCC':'A', 'GCA':'A', 
                         'GCG':'A', 'GAU':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 
                         'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'};
  var aa = [];

  for (i = 0; i < rna.length; i+= 3) {
    aa.push(rna_codon_table[rna.substring(i, i+3)]);
  }
  
  return aa.join('');
}

function kmers(k, dna) {
  
  var list = [];
  for (i = 0; i < dna.length - k + 1; i++) {
    list.push(dna.substring(i, i + k));
  }
  return list;
}

function DeBruijn(k, dna) {
    
  var nodes = [];
  var ins = [];
  var outs = [];
  var temp = [];
  nodes = kmers(k, dna);
  outs = kmers(k-1, dna.substring(1, dna.length));
  ins = kmers(k-1, dna.substring(0, dna.length - 1));
  
  nodes.sort();
  ins.sort();
  outs.sort();
   
  for (i = 0; i < ins.length; i++) {
    for (j = 0; j < outs.length; j++) {
      var node = ins[i].concat(outs[j].substring(k - 2, k));
      var index = nodes.indexOf(node);
      if (index != -1){
        var arr = [ins[i], node, outs[j], '<br>'];
        temp.push(arr);
        nodes.splice(index, 1);
      }
    }
  }
  
  return temp;
}

function linearSpectrum(peptide) {
  var amino_acid_mass = { 'G' : 57,  'A' : 71,  'S' : 87,  'P' : 97, 'V' : 99,  'T' : 101, 'C' : 103, 'I' : 113, 'L' : 113, 'N' : 114, 'D' : 115, 'K' : 128, 'Q' : 128, 'E' : 129,'M' : 131, 'H' : 137, 'F' : 147, 'R' : 156, 'Y' : 163, 'W' : 186};
  var prefix_mass = [0];
  
  for (i = 1; i <= peptide.length; i++){
    prefix_mass.push(prefix_mass[i-1] + amino_acid_mass[peptide[i-1]]);
  }
  
  for (i = 1; i <= peptide.length; i++){
    for (j = i + 1; j <= peptide.length; j++){
      prefix_mass.push(prefix_mass[j] - prefix_mass[i]);
    }
  }
  
  prefix_mass.sort(function(a, b) { return a - b; });
  
  return prefix_mass;
}

var a = DeBruijn(number, text);

myHeading.innerHTML = a;
