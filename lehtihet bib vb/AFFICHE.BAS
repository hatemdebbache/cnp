Attribute VB_Name = "AFFICHE"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << AFFICHE >>
'  Requis   :
'  Externes :
'
'  Descrip. : Affichages de tableaux
'
Option Explicit
Sub Vecteur(ByVal k1 As Integer, _
            ByVal k2 As Integer, _
            ByRef VEC() As Double, _
   Optional ByVal f As String = "")
    
    Dim i As Integer, k As Integer
    Dim N As Integer, s As Integer
    Dim v As Double, ff As String
    
    s = Sgn(k2 - k1)      ' Sens d'affichage
    N = s * (k2 - k1) + 1 ' Nbr. d'éléments à afficher
    k = k1                ' Premier élément
    For i = 1 To N
        v = VEC(k)
        If v >= 0 Then ff = " " Else ff = ""
        Debug.Print ff + Format(v, f)
        k = k + s         ' Avance/Recule selon S
    Next i
End Sub
Sub Matrice(ByVal ik1 As Integer, _
            ByVal ik2 As Integer, _
            ByVal jk1 As Integer, _
            ByVal jk2 As Integer, _
            ByRef MAT() As Double, _
   Optional ByVal f As String = "")
    
    Dim i As Integer, ik As Integer
    Dim j As Integer, jk As Integer
    Dim Si As Integer, Ni As Integer
    Dim Sj As Integer, Nj As Integer
    Dim v As Double, ff As String
    
    Si = Sgn(ik2 - ik1): Ni = Si * (ik2 - ik1) + 1
    Sj = Sgn(jk2 - jk1): Nj = Sj * (jk2 - jk1) + 1
    ik = ik1
    For i = 1 To Ni
        jk = jk1
        For j = 1 To Nj
            v = MAT(ik, jk)
            If v >= 0 Then ff = " " Else ff = ""
            Debug.Print ff + Format(v, f);
            jk = jk + Sj
        Next j
        Debug.Print
        ik = ik + Si
    Next i
End Sub



