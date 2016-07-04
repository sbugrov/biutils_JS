var myHeading = document.getElementById("demo");

var text = 'TAATGCCATGGGATGTT';
var number = 3;


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

var a = DeBruijn(number, text);

myHeading.innerHTML = a;

// |,AA,AAT,AT,|,
// |,AT,ATG,GG,|,
// |,AT,ATG,GG,|,
// |,AT,ATG,TG,|,
// |,CA,CAT,AT,|,
// |,CC,CCA,AA,|,
// |,GA,GAT,AT,|,
// |,GC,GCC,CC,|,
// |,GG,GGA,AA,|,
// |,GG,GGG,GG,|,
// |,GT,GTT,AT,|,
// |,TA,TAA,AA,|,
// |,TG,TGT,AT,|,
// |,TG,TGC,CC,|,
// |,TG,TGG,GG,|

