/*
 * The MIT License
 *
 * Copyright 2016 Mikl&oacute;s Cs&#369;r&ouml;s.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package splice.model;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class TreeNode 
{
    private String      name;
    private TreeNode    mamma;          // parent node
    private int         num_children;   // number of children
    private TreeNode[]  child;          // the children: left to right
    private double      length;         // edge length to parent
    private int child_capacity; //

    private int id;
    private static int next_id=1;
    public TreeNode()
    {
        child_capacity=num_children=0;
        setId(next_id++);
    }
  
    public final int getId(){return id;}
    public final void setId(int id){this.id=id;}

    public String getName(){return name;}
    public void setName(String s){name=s;}

    public double getLength(){return length;}
    public void setLength(double l){length=l;}

    public TreeNode getParent(){return mamma;}
    protected void setParent(TreeNode n){mamma=n;}

    public int getNumChildren(){return num_children;}

    /**
     * @param idx
     * @return null if no child with that index
     */
    public TreeNode getChild(int idx)
    {
        if (idx>=0 && idx<num_children)
            return child[idx];
        else return null;
    }

    public static final int LEFT_CHILD=0;
    public static final int RIGHT_CHILD=1;
  
    /**
     * Sets the child at the given index to n: also sets the parent of n to be this.
     * @param n
     * @param child_index
     */
    public void setChild(TreeNode n, int child_index)
    {
        child[child_index]=n;
        n.mamma=this;
        n.index_at_parent=child_index;
    }

    public void clearChildren()
    {
        num_children=0;
    }
  
    private int index_at_parent=-1;
    /**
     * Adds a new child; the child's parent is set to this node
     * @param n
     */
    public void addChild(TreeNode n)
    {
        // expand array if necessary
        if (num_children==child_capacity)
        {
            int new_capacity=(child_capacity==0?2:child_capacity*2);
            TreeNode[] new_child=new TreeNode[new_capacity];
            if (child_capacity != 0)
              System.arraycopy(child,0,new_child,0,child_capacity);
            child=new_child;
            child_capacity=new_capacity;
        }
        setChild(n, num_children);
        num_children++;
    }
  
    public TreeNode newChild()
    {
        TreeNode n=new TreeNode();
        this.addChild(n);
        return n;
    }
    public boolean isRoot(){return mamma==null;}
    public boolean isLeaf(){return num_children==0;}

    public int getIndexAtParent(){return index_at_parent;}
    public boolean isOnTheLeft(){return index_at_parent==LEFT_CHILD;}
    public boolean isOnTheRight(){return index_at_parent==RIGHT_CHILD;}

    public static boolean ALWAYS_QUOTE_NAME=false;
    private static final String EMPTY_QUOTED_NAME = ""+Parser.QUOTE+Parser.QUOTE;

    public String newickName(){return newickName(ALWAYS_QUOTE_NAME);}
    public String newickName(boolean always_quote)
    {
        String nnn = getName();
        if (nnn==null) return (always_quote?EMPTY_QUOTED_NAME:"");

        char[] cname=nnn.toCharArray();
        boolean need_quote=always_quote;
        if (!always_quote)
          for (int i=0; i<cname.length; i++)
            if (Parser.NEED_QUOTE_FOR.indexOf(cname[i])!=-1 || Character.isWhitespace(cname[i]))
            {
              need_quote=true;
              break;
            }
        if (need_quote)
        {
          StringBuilder sb=new StringBuilder();

          sb.append(Parser.QUOTE);
          int from_idx=0;
          int to_idx;
          do {
            to_idx=name.indexOf(Parser.QUOTE,from_idx);
            if (to_idx==-1)
              sb.append(cname,from_idx,cname.length-from_idx);
            else
            {
              sb.append(cname,from_idx,to_idx-from_idx+1);
              sb.append(Parser.QUOTE);
              from_idx=to_idx+1;
            }
          } while(to_idx != -1 && from_idx<cname.length);
          sb.append(Parser.QUOTE);
          return sb.toString();
        } else
          return nnn;
    }

    public static boolean SHOW_EDGE_LENGTHS=true;
    public static boolean LINE_BREAKS_AFTER_EACH_NODE=false;
    public static boolean SHOW_NODE_ID_INFO=false;

    public String newickTree()
    {
      return newickTree(ALWAYS_QUOTE_NAME,SHOW_EDGE_LENGTHS,LINE_BREAKS_AFTER_EACH_NODE,SHOW_NODE_ID_INFO);
    }

    public String newickTree(
      boolean always_quote_name,
      boolean show_edge_lengths,
      boolean line_breaks_after_each_node,
      boolean show_node_id_info)
    {
      return newickSubtree(always_quote_name,show_edge_lengths,line_breaks_after_each_node,show_node_id_info)+";";
    }

    public String newickSubtree()
    {
      return newickSubtree(ALWAYS_QUOTE_NAME,SHOW_EDGE_LENGTHS,LINE_BREAKS_AFTER_EACH_NODE,SHOW_NODE_ID_INFO);
    }

    public String newickSubtree(
      boolean always_quote_name,
      boolean show_edge_lengths,
      boolean line_breaks_after_each_node,
      boolean show_node_id_info)
    {
      //System.out.println("#**TN.nS "+toString());

      StringBuilder sb=new StringBuilder();
      if (show_node_id_info)
      {
        sb.append(Parser.LBRACKET);
        sb.append("<<< ");
        sb.append(id);
        sb.append(Parser.RBRACKET);
      }
      if (num_children>0)
      {
        sb.append(Parser.LPAREN);
        for (int i=0; i<num_children; i++)
        {
          if (i>0)
          {
            sb.append(Parser.COMMA);
            if (line_breaks_after_each_node)
              sb.append("\n");
            else
              sb.append(' ');
            if (show_node_id_info)
            {
              sb.append(Parser.LBRACKET);
              sb.append("=== ");
              sb.append(id);
              sb.append(Parser.RBRACKET);
            }
          }
          sb.append(child[i].newickSubtree(always_quote_name,show_edge_lengths,line_breaks_after_each_node,show_node_id_info));
        }
        sb.append(Parser.RPAREN);
      }

      if (isLeaf() || getName() != null)
      {
          if (!isLeaf())
              sb.append(' ');
          sb.append(newickName(always_quote_name));
      }

      if (show_edge_lengths && !isRoot())
      {
        sb.append(Parser.COLON);
        sb.append(newickEdgeLength());
      }

      if (show_node_id_info)
      {
        sb.append(Parser.LBRACKET);
        sb.append(" >>> ");
        sb.append(id);
        sb.append(Parser.RBRACKET);
      }


      return sb.toString();
    }
  
    protected String newickEdgeLength()
    {
        return toIEEEFormat(length);
    }

    public String shortDesc()
    {
        String nnn = getName();
        return "$"+id+"_"+(nnn==null?"noname":nnn);
    }

    public String subtreeParamString()
    {
      if (isLeaf())
        return toString();
      StringBuilder sb=new StringBuilder(toString());
      for (int i=0; i<num_children; i++)
      {
        sb.append('\n');
        sb.append(child[i].subtreeParamString());
      }
      return sb.toString();
    }
  
    protected String paramString()
    {
      StringBuilder sb=new StringBuilder();
      sb.append(shortDesc());
      sb.append(" len ");
      sb.append(length);
      sb.append(" prnt ");
      if (mamma != null)
        sb.append(mamma.shortDesc());
      else
        sb.append('-');
      sb.append(" chld [");
      for (int i=0; i<num_children; i++)
      {
        if (i>0)
          sb.append(", ");
        sb.append(child[i].shortDesc());
      }
      sb.append("]");
      return sb.toString();
    }

    private String class_id_name = null;
    @Override
    public String toString()
    {

      if (class_id_name == null)
      {
          String simple_name = getClass().getSimpleName();
          class_id_name = "";
          for (int i=0; i<simple_name.length(); i++)
          {
              char c = simple_name.charAt(i);
              if (Character.isUpperCase(c))
                  class_id_name = class_id_name + c;
          }
          if ("".equals(class_id_name))
              class_id_name = simple_name.substring(0,Math.min(simple_name.length(),5));
      }
      StringBuilder sb=new StringBuilder(class_id_name);
      sb.append("[");
      sb.append(paramString());
      sb.append("]");
      return sb.toString();
    }

    static int DECIMALS=4;
    private static final double rounding_factor = 1000.;
    private static final double LOG10=Math.log(10.);
    private static final double short_branch = 1e-7;

    private static final double[] EXP10 = {1.,10.,100.,1000.,10000.,100000.,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13};

    public static String toIEEEFormat(double d)
    {

      if (Double.isNaN(d))
        return "NaN";
      else if (d==Double.POSITIVE_INFINITY)
        return "Inf";
      else if (d==Double.NEGATIVE_INFINITY)
        return "-Inf";
      else if (d==0.)
        return "0";
      else
      {
        int magnitude=(int) (Math.log(d)/LOG10+0.01);
        String retval;
        if (magnitude>=-DECIMALS)
        {
          double r=((int)(d*rounding_factor+0.5))/rounding_factor;
          retval= r+"";
        } else if (d<short_branch) {
            retval= "9.9e-9";//(short_branch*0.5)+"";
        } else
        {
            // only one digit after decimal point
            double m = EXP10[-magnitude-1];
            double r=((int)(d*m)+0.5)/m;

            retval= r+"";
        }
        //System.out.println("#**TN.tIEEE "+d+"\tmag "+magnitude+"\tr "+retval);
        return retval;
        //return d+"";
      }
    }

//    /**
//     * Calculates an array of nodes in the subtree rooted at this node.
//     * The nodes are listed in postorder traversal.
//     *
//     * (The node itself is always included.)
//     */
//    public TreeNode[] subtreeNodes()
//    {
//        List<TreeNode> nodes_in_subtree = new ArrayList<>();
//        subtreeNodes(this,nodes_in_subtree,false);
//
//        return nodes_in_subtree.toArray(new TreeNode[0]);
//    }
//  
//    /**
//     * Calculates an array of leaves in the subtree rooted at this node.
//     * (The node itself is included if it happens to be a leaf.)
//     */
//    public TreeNode[] subtreeLeaves()
//    {
//        List<TreeNode> leaves_in_subtree = new ArrayList<>();
//        subtreeNodes(this,leaves_in_subtree,true);
//
//        return leaves_in_subtree.toArray(new TreeNode[0]);
//    }
  
    /**
     * Auxiliary recursive function for collecting nodes in the subtree by performing a 
     * postorder traversal.
     *
     * @param n subtree root
     * @param nodes_so_far vector to which the nodes are added in the order they are visited
     * @param leaves_only whether only leaves should be collected or all nodes
     */
    private void subtreeNodes(TreeNode n, List<TreeNode> nodes_so_far, boolean leaves_only)
    {
        if (n.isLeaf())
            nodes_so_far.add(n);
        else
        {
            for (int cidx=0; cidx<n.getNumChildren(); cidx++)
               subtreeNodes(n.getChild(cidx),nodes_so_far,leaves_only);
            if (!leaves_only) nodes_so_far.add(n);
        }

    }

//    // used for retrieving the path to the root
//    private Stack<TreeNode> stack=new Stack<TreeNode>();
//    private void getPathToRoot(Stack<TreeNode> so_far)
//    {
//      so_far.push(this);
//      if (mamma != null)
//        mamma.getPathToRoot(so_far);
//    }
//
//    /** 
//     * @return a Vector of TreeNodes with current object in position 0, then its parent in pos 1, 
//     * and so on until the root in the last position
//     */
//    public List<TreeNode> getPathToRoot()
//    {
//      stack.clear();
//      getPathToRoot(stack);
//      List<TreeNode> stack_copy = new ArrayList<>(stack.size());
//      stack_copy.addAll(stack);
//      return stack_copy;
//    }
//
//    /**
//     * Computes the lowest common ancestor between this guy and another node.
//     * There is no special data structure, simply the paths from the root are matched
//     * until the first discrepancy is found.
//     */
//    public TreeNode getLowestCommonAncestorWith(TreeNode n)
//    {
//      this.stack.clear();
//      n.stack.clear();
//      this.getPathToRoot(this.stack);
//      n.getPathToRoot(n.stack);
//
//      TreeNode lca = (TreeNode)this.stack.pop(); // root: common ancestor by def
//      n.stack.pop();
//
//      while (!this.stack.empty() && !n.stack.empty())
//      {
//        TreeNode on_path=(TreeNode)this.stack.pop();
//        if (on_path.equals(n.stack.pop()))
//          lca=on_path;
//        else
//          break;
//      }
//
//      return lca;
//    }
//  
    /** 
     * Re-roots the tree at this node
     */
    public void reroot()
    {
        if (isRoot())
            return;
        TreeNode p = getParent();
        p.reroot();
        int idx = getIndexAtParent();
        for (int i=idx; i+1<p.getNumChildren(); i++)
            p.setChild(p.getChild(i+1), i);
        p.num_children--;
        addChild(p);
        p.setLength(getLength());
        setLength(0.);
        setParent(null);
    }
  
    /**
     * Number of nodes in the subtree rooted here (including this very node)
     */
    public int numNodes()
    {
        int s = 0;
        for (int cidx=0; cidx<getNumChildren(); cidx++)
            s+=getChild(cidx).numNodes();
        return s+1;
    }

    @Override
    public boolean equals(Object o)
    {
      if (o==null)
        return false;
      return id==((TreeNode)o).id;
    }

    /**
     * Gives back a Traversal object for the subtree rooted here.
     */
    public Traversal getTraversal()
    {
        return new Traversal(this);
    }

//    public static void main(String[] args) throws Exception
//    {
//        if (args.length <1)
//        {
//            throw new IllegalArgumentException("Call as "+TreeNode.class.getCanonicalName()+" treefile");
//        }
//
//        if ("epslevels".equals(args[0]))
//        {
//            TreeNode root = Parser.readNewick(new java.io.BufferedReader(new java.io.FileReader(args[1])));
//            TreeNode[] nodes = root.getTraversal().getDFT();
//
//            Map<TreeNode,Integer> node_indexes = new HashMap<TreeNode,Integer>();
//            Map<TreeNode, Integer> node_heights = new HashMap<TreeNode, Integer>();
//
//            for (int node_idx=0; node_idx<nodes.length; node_idx++)
//            {
//                TreeNode N = nodes[node_idx];
//                node_indexes.put(N, node_idx);
//                if (N.isLeaf())
//                {
//                    node_heights.put(N, 0);
//                } else
//                {
//                    int nc = N.getNumChildren();
//                    int max_hgt = 0;
//                    for (int ci=0; ci<nc; ci++)
//                    {
//                        TreeNode C = N.getChild(ci);
//                        int h = node_heights.get(C);
//                        if (h>max_hgt) max_hgt=h;
//                    }
//                    node_heights.put(N, (max_hgt+1));
//                }
//            }
//            System.out.println("% ca.umontreal.iro.evolution.TreeNode "+args[0]+"\t"+args[1]);
//
//            System.out.println("[ %NAME");
//            for (int node_idx=0; node_idx<nodes.length; node_idx++)
//                System.out.println("   ("+nodes[node_idx].newickName()+")\t % "+node_idx);
//            System.out.println("] /node.names exch def");
//
//            System.out.println("[ %HEIGHT");
//            for (int node_idx=0; node_idx<nodes.length; node_idx++)
//                System.out.println("   "+node_heights.get(nodes[node_idx])+"\t % "+node_idx+"/"+nodes[node_idx].newickName());
//            System.out.println("] /node.heights exch def");
//
//            System.out.println("[ %PARENT");
//            for (int node_idx=0; node_idx<nodes.length; node_idx++)
//            {
//                TreeNode N = nodes[node_idx];
//                int pidx = (N.isRoot()?-1:node_indexes.get(N.getParent()));
//                System.out.println("   "+pidx+"\t % "+node_idx+"/"+nodes[node_idx].newickName());
//            }
//            System.out.println("] /node.parents exch def");
//
//        } else
//        {
//            TreeNode root = Parser.readNewick(new java.io.BufferedReader(new java.io.FileReader(args[0])));
//            TreeNode[] nodes = root.getTraversal().getDFT();
//            Map<TreeNode, Double> position = new HashMap<>();
//            for (int node_idx=nodes.length-1; node_idx>=0; node_idx--)
//            {
//                TreeNode N = nodes[node_idx];
//                if (N.isRoot())
//                    position.put(N, 0.0);
//                else
//                {
//                    TreeNode P = N.getParent();
//                    double Ppos = position.get(P).doubleValue();
//                    double Npos = Ppos + 450 * N.getLength();
//                    position.put(N, Npos);
//                }
//            }
//            System.out.println("Node\tcumulative depth\tlength");
//            for (int node_idx=0; node_idx<nodes.length; node_idx++)
//            {
//                TreeNode N = nodes[node_idx];
//                System.out.println(node_idx+"/"+N.newickName()+"\t"+position.get(N)+"\t"+N.getLength());
//            }
//        }
//    }
  

    public class Traversal 
    {
      public Traversal(TreeNode subtree_root)
      {
          this.subtree_root = subtree_root;
      }

      private TreeNode subtree_root;

      /** 
       * Puts the nodes in an array suitable for depth-first traversal of the subtree rooted dft_root.
       */
      public TreeNode[] getDFT()
      {
          List<TreeNode> dft_vector=new ArrayList<>();
          collectNodes(subtree_root, dft_vector, false);
          return dft_vector.toArray(new TreeNode[0]);
      }

      /** 
       * @return array of leaves of the subtree rooted at this guy.
       */
      public TreeNode[] getLeaves()
      {
          List<TreeNode> leaf_vector = new ArrayList<>();
          collectNodes(subtree_root, leaf_vector, true);
          return leaf_vector.toArray(new TreeNode[0]);
      }

      /**
       * Helper recursive function for filling up the vector of nodes for depth-first traversal
       */
      private void collectNodes(TreeNode N, List<TreeNode> V, boolean leaf_only)
      {
          for (int child_idx=0; child_idx<N.getNumChildren(); child_idx++){
              TreeNode Child = (TreeNode)N.getChild(child_idx);
              collectNodes(Child,V,leaf_only);
          }
          if (!leaf_only || N.isLeaf())
              V.add(N);
      }
    }
}
